# KIT-012: what the DA step leaves resident in the application process at 300 samples.
# The KIT-011 candidate's trace (kit-journey-011/candidate-1, the on-disk default) ends at 3.02 GiB anonymous
# memory after the DA step against 2.26 GiB after the session export, so the DA step accounts for most of
# the resident memory above the pre-DA level. This probe replays the DA tab's production path in process on
# that member's project (reflink copy): session load through the native Open, the reference DA analysis, the
# results-table page, the static volcano, the interactive (Glimma) volcano and the heatmap, as the journey
# driver exercises them. At every boundary it records the process's anonymous resident memory as sampled
# (raw), after a full collection and after the allocator released its free pages, together with the live R
# heap and the size and kind (reference or payload) of every module holder, then releases every holder and
# measures again. The split tells the ticket whether a holder or allocator retention owns that memory.
# Not a counted run (no stage, no ledger line); the record is written beside this script.
# Usage: R_LIBS=<proof>/core/browser-smoke.X9aj9P/library Rscript docs/qa/workflow-execution/kit/kit-012/da-resident-probe.R
# With KIT012_PROBE_WALK=1 the probe also walks every loaded namespace, the store and the session for live
# bytes (da-resident-probe-walk.json). The walk forces the lazy-loaded data of namespaces the DA step loads,
# which adds about 0.35 GiB of live heap of its own, so the clean record is the run without it.
suppressMessages(devtools::load_all(quiet = TRUE))
repo <- normalizePath(".")
h <- new.env(parent = globalenv())
sys.source(file.path(repo, "tools/profiling/workflow_kit.R"), h)
b <- h$kitBindings(file.path(repo, "docs/qa/workflow-execution/kit/bindings.json"))
member <- file.path(b$root, "kit-journey-011", "candidate-1")
journey <- readRDS(file.path(member, "journey-result.rds"))
stopifnot(identical(journey$status, "complete"))
work <- file.path(Sys.getenv("HOME"), ".cache", "multischolar-kit012-da-probe")
unlink(work, recursive = TRUE)
dir.create(work, recursive = TRUE)
status <- system2("cp", c("-a", "--reflink=auto", file.path(member, "project", "dia_baseline_smoke"), work))
stopifnot(status == 0L)
on.exit(unlink(work, recursive = TRUE), add = TRUE)
project <- file.path(work, "dia_baseline_smoke")
results <- file.path(project, "results", "proteomics")
paths <- list(base_dir = project, results_dir = results, source_dir = file.path(project, "scripts", "proteomics"),
    protein_qc_dir = file.path(results, "protein_qc"), peptide_qc_dir = file.path(results, "peptide_qc"),
    da_output_dir = file.path(results, "da_proteins"), publication_graphs_dir = file.path(results, "publication_graphs"),
    time_dir = file.path(results, "publication_graphs"), integration_dir = file.path(project, "integration"))
for (path in unlist(paths)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
source_dir <- normalizePath(paths$source_dir, mustWork = TRUE)
archive <- list.files(source_dir, "[.]msdia[.]tar$", full.names = TRUE)
stopifnot(length(archive) == 1L)

gib <- function(bytes) round(bytes / 2^30, 3)
live_heap <- function(collected) sum(collected[c("Ncells", "Vcells"), "used"] * c(56, 8))
anon <- function() snapshotWorkflowMemory()$process$RssAnon
holder_kind <- function(value) {
    if (is.null(value)) return("NULL")
    if (isSavedWorkflowReference(value)) return("reference")
    if (isS4(value)) return(paste("S4", class(value)))
    if (is.data.frame(value)) return(sprintf("data.frame %d x %d", nrow(value), ncol(value)))
    class(value)[[1L]]
}
describe <- function(owner, fields) {
    stats::setNames(lapply(fields, function(field) {
        value <- shiny::isolate(owner[[field]])
        list(kind = holder_kind(value), bytes = as.numeric(utils::object.size(value)))
    }), fields)
}
# The largest objects an environment holds (what the live R heap is made of).
inventory <- function(env, limit = 12L) {
    names <- ls(env, all.names = TRUE)
    if (!length(names)) return(list())
    bytes <- vapply(names, function(name) {
        value <- tryCatch(get(name, envir = env, inherits = FALSE), error = function(e) NULL)
        if (is.environment(value)) return(NA_real_)
        as.numeric(utils::object.size(value))
    }, numeric(1))
    keep <- order(bytes, decreasing = TRUE, na.last = TRUE)[seq_len(min(limit, length(names)))]
    stats::setNames(as.list(bytes[keep]), names[keep])
}
# Bytes of the non-environment bindings an environment holds.
env_bytes <- function(env) {
    names <- ls(env, all.names = TRUE)
    if (!length(names)) return(0)
    sum(vapply(names, function(name) {
        value <- tryCatch(get(name, envir = env, inherits = FALSE), error = function(e) NULL)
        if (is.environment(value)) 0 else as.numeric(utils::object.size(value))
    }, numeric(1)))
}
# Where live R memory can hide beyond the module holders: the package's namespace-level environments, every
# other loaded namespace's top-level bindings, the entity store's record and payload caches, the session's
# user data, Shiny's globals and the options.
reachable <- function(store, session) {
    out <- list()
    ns <- asNamespace("MultiScholaR")
    names <- ls(ns, all.names = TRUE)
    for (name in names[!startsWith(names, ".__")]) {
        value <- get(name, envir = ns, inherits = FALSE)
        if (is.environment(value)) out[[paste0("MultiScholaR/", name)]] <- env_bytes(value)
    }
    for (name in setdiff(loadedNamespaces(), "MultiScholaR")) {
        bytes <- tryCatch(env_bytes(asNamespace(name)), error = function(e) NA_real_)
        if (is.finite(bytes) && bytes > 16 * 2^20) out[[paste0("namespace/", name)]] <- bytes
    }
    out[["store/cache_bytes"]] <- store$cache_bytes
    out[["store/records"]] <- env_bytes(store$records)
    out[["store/cache"]] <- env_bytes(store$cache)
    out[["store/fields"]] <- env_bytes(store)
    out[["session/userData"]] <- env_bytes(session$userData)
    out[["shiny/.globals"]] <- env_bytes(get(".globals", envir = asNamespace("shiny")))
    options_bytes <- vapply(options(), function(value) as.numeric(utils::object.size(value)), numeric(1))
    out[["options/largest"]] <- as.list(sort(options_bytes, decreasing = TRUE)[1:3])
    collected <- gc(full = TRUE)
    out[["gc/ncells_used"]] <- collected["Ncells", "used"]
    out[["gc/vcells_used"]] <- collected["Vcells", "used"]
    out
}
walk <- identical(Sys.getenv("KIT012_PROBE_WALK"), "1")
record_path <- file.path(repo, "docs/qa/workflow-execution/kit/kit-012",
    if (walk) "da-resident-probe-walk.json" else "da-resident-probe.json")
record <- list(schema = "kit.da_resident_probe.v1", mode = if (walk) "walk" else "clean", ticket = "KIT-012", date = format(Sys.Date()),
    head_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
    member = "kit-journey-011/candidate-1 (kit_011, 300-sample corpus, on-disk default)",
    trace = list(anon_after_session_export_gib = 2.26, anon_after_da_gib = 3.09, resident_end_gib = 3.017),
    allocator_release_supported = isTRUE(workflowAllocatorRelease(FALSE)), boundaries = list())
boundary <- function(label, seconds = NA_real_, holders = NULL) {
    raw <- anon()
    collected <- gc(full = TRUE)
    after_gc <- anon()
    released <- workflowAllocatorRelease(TRUE)
    after_release <- anon()
    entry <- list(label = label, seconds = seconds, anon_raw = raw, anon_after_gc = after_gc,
        anon_after_release = after_release, live_r_heap = live_heap(collected), pages_released = released,
        holders = holders)
    record$boundaries[[length(record$boundaries) + 1L]] <<- entry
    cat(sprintf("%-22s raw=%.3f gc=%.3f trim=%.3f liveR=%.3f GiB%s\n", label, gib(raw), gib(after_gc),
        gib(after_release), gib(entry$live_r_heap), if (is.na(seconds)) "" else sprintf("  (%.1f s)", seconds)))
    invisible(entry)
}
timed <- function(expr) {
    started <- proc.time()[["elapsed"]]
    value <- force(expr)
    list(value = value, seconds = round(proc.time()[["elapsed"]] - started, 1))
}

boundary("start")
roots <- list.files(source_dir, "^[.]dia-native-[0-9a-f]+$", all.files = TRUE, include.dirs = TRUE, full.names = TRUE)
stopifnot(length(roots) == 1L)
identity <- readSavedStoreIdentity(file.path(roots, ".workflow-entities", "store.json"))
manager <- newDIAStateAdapter(newWorkflowEntityStore(roots, identity$project_id, identity$workflow_id,
    max_cache_bytes = 0), reopen = TRUE)
manager$.context$session_format <- "native_v1"
manager$.context$reference_states <- TRUE
manager$.context$experiment_paths <- paths
workflow <- shiny::reactiveValues(state_manager = manager, data_type = "peptide", data_format = "diann",
    tab_status = list(), contrasts_tbl = NULL, uniprot_dat_cln = NULL, design_matrix = NULL, config_list = NULL,
    da_analysis_results_list = NULL, da_ui_params = NULL)
bus_fields <- c("contrasts_tbl", "uniprot_dat_cln", "design_matrix", "config_list", "da_analysis_results_list")
da_fields <- c("current_s4_object", "da_results_list", "current_heatmap_plot", "current_row_clusters", "current_col_clusters")
probe_error <- NULL
tryCatch(shiny::testServer(function(input, output, session) {}, {
    workflow$state_manager <- attachDIAStateSession(manager, session, workflow)
    boundary("store_open")
    # The DA tab's "Load filtered session" in a fresh process: the native Open of the journey's export and its
    # adoption as the session owner (the journey itself resolved its own export in-session).
    loaded <- timed({
        imported <- openDIANativeSession(archive, source_dir)
        adoptDIANativeSession(imported, workflow, session, source_dir, paths)
    })
    view <- loaded$value
    stopifnot(isSavedWorkflowReference(view$current_s4_object))
    workflow$contrasts_tbl <- view$contrasts_tbl
    workflow$design_matrix <- view$design_matrix
    workflow$config_list <- view$config_list
    uniprot_path <- file.path(source_dir, "uniprot_dat_cln.RDS")
    workflow$uniprot_dat_cln <- if (is.data.frame(view$uniprot_dat_cln)) view$uniprot_dat_cln else
        if (file.exists(uniprot_path)) readRDS(uniprot_path) else NULL
    da_data <- shiny::reactiveValues(current_s4_object = view$current_s4_object,
        contrasts_available = view$contrasts_tbl$contrasts, analysis_complete = FALSE, da_results_list = NULL,
        current_heatmap_plot = NULL, current_row_clusters = NULL, current_col_clusters = NULL)
    view_fields <- names(view)
    view_sizes <- vapply(view_fields, function(f) as.numeric(utils::object.size(view[[f]])), numeric(1))
    record$session_view <<- list(fields = view_fields, bytes = unname(view_sizes),
        kinds = unname(vapply(view_fields, function(f) holder_kind(view[[f]]), character(1))),
        current_state = workflow$state_manager$current_state, states = length(workflow$state_manager$.context$snapshot$states))
    view <- imported <- NULL
    boundary("session_loaded", loaded$seconds, list(bus = describe(workflow, bus_fields), da = describe(da_data, da_fields)))
    if (walk) record$reachable_after_session_load <<- reachable(workflow$state_manager$.context$store, session)
    settings <- list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0)
    analysed <- timed(runDIAReferenceDaAnalysis(settings, da_data, workflow, paths))
    manifest <- analysed$value
    stopifnot(isDIAResultManifest(manifest), identical(shiny::isolate(da_data$da_results_list), manifest))
    record$da <<- list(contrasts = manifest$contrasts$friendly_names, results = length(manifest$results),
        manifest_bytes = as.numeric(utils::object.size(manifest)), seconds = analysed$seconds)
    contrast <- manifest$contrasts$friendly_names[[1L]]
    boundary("da_analysis", analysed$seconds, list(bus = describe(workflow, bus_fields), da = describe(da_data, da_fields)))
    record$inventory_after_da <<- list(adapter_context = inventory(workflow$state_manager$.context),
        probe_frame = inventory(environment()), global = inventory(globalenv()))
    if (walk) record$reachable_after_da <<- reachable(workflow$state_manager$.context$store, session)
    page <- timed(readDIAResultPage(manifest, contrast, "all", 0.05, 0, 0L, 100L))
    record$table_page <<- list(rows = nrow(page$value$data), full_population_rows = page$value$full_population_rows,
        seconds = page$seconds)
    page <- NULL
    boundary("results_table_page", record$table_page$seconds)
    static <- timed({
        plot <- generateProtDAVolcanoStatic(manifest, contrast, da_q_val_thresh = 0.05, lfc_threshold = 0,
            show_labels = TRUE, n_labels = 10)
        device <- file.path(work, "static-volcano.png")
        grDevices::png(device, width = 1200, height = 900)
        print(plot)
        grDevices::dev.off()
        as.numeric(utils::object.size(plot))
    })
    record$static_volcano <<- list(plot_bytes = static$value, seconds = static$seconds)
    boundary("static_volcano", static$seconds)
    glimma <- timed({
        widget <- generateProtDAVolcanoPlotGlimma(manifest, contrast, da_q_val_thresh = 0.05,
            args_row_id = diaGuiProteinKey(shiny::isolate(da_data$current_s4_object)), output_dir = paths$da_output_dir)
        stopifnot(!is.null(widget))
        as.numeric(utils::object.size(widget))
    })
    record$interactive_volcano <<- list(widget_bytes = glimma$value, seconds = glimma$seconds)
    boundary("interactive_volcano", glimma$seconds)
    heat <- timed({
        result <- generateProtDAHeatmap(manifest, contrast, top_n_genes = 50, clustering_method = "ward.D2",
            distance_method = "euclidean", cluster_rows = TRUE, cluster_cols = TRUE, scale_data = "row",
            color_scheme = "RdBu", show_gene_names = FALSE, da_q_val_thresh = 0.05, tree_cut_method = "none",
            n_clusters = 4, cut_height = 0.5, min_cluster_size = 3)
        # The module keeps the plot and the cluster vectors on its reactive owner; with no protein under the
        # q-value threshold the heatmap is NULL (the synthetic corpus has no differential proteins).
        if (is.list(result) && !is.null(result$plot)) {
            da_data$current_row_clusters <- result$row_clusters
            da_data$current_col_clusters <- result$col_clusters
            da_data$current_heatmap_plot <- result$plot
            device <- file.path(work, "heatmap.png")
            grDevices::png(device, width = 1200, height = 900)
            print(result$plot)
            grDevices::dev.off()
            as.numeric(utils::object.size(result$plot))
        } else NA_real_
    })
    record$heatmap <<- list(plot_bytes = heat$value, seconds = heat$seconds,
        note = if (is.na(heat$value)) "no protein under the q-value threshold: the heatmap is NULL" else NULL)
    boundary("heatmap", heat$seconds, list(bus = describe(workflow, bus_fields), da = describe(da_data, da_fields)))
    # Every module holder released: what remains is not a holder.
    for (field in da_fields) da_data[[field]] <- NULL
    for (field in bus_fields) workflow[[field]] <- NULL
    manifest <- NULL
    boundary("holders_released", holders = list(bus = describe(workflow, bus_fields), da = describe(da_data, da_fields)))
    record$inventory_after_release <<- list(adapter_context = inventory(workflow$state_manager$.context),
        probe_frame = inventory(environment()), global = inventory(globalenv()))
    if (walk) {
    record$reachable_after_release <<- reachable(workflow$state_manager$.context$store, session)
    changed <- Filter(function(name) !identical(record$reachable_after_release[[name]], record$reachable_after_session_load[[name]]),
        names(record$reachable_after_release))
    for (name in changed) cat(sprintf("reachable %-40s %s -> %s\n", name,
        paste(format(unlist(record$reachable_after_session_load[[name]]), big.mark = ","), collapse = "/"),
        paste(format(unlist(record$reachable_after_release[[name]]), big.mark = ","), collapse = "/")))
    }
    closeDIAStateAdapter(workflow$state_manager)
    boundary("store_closed")
}), error = function(error) { probe_error <<- conditionMessage(error); cat("PROBE ERROR:", probe_error, "\n") })
record$error <- probe_error
first <- function(label) {
    found <- Filter(function(e) identical(e$label, label), record$boundaries)
    if (!length(found)) return(list(anon_raw = NA_real_, anon_after_release = NA_real_, live_r_heap = NA_real_))
    found[[1L]]
}
delta <- function(from, to, field) first(to)[[field]] - first(from)[[field]]
record$result <- list(
    da_step_anon_raw_delta = delta("session_loaded", "heatmap", "anon_raw"),
    da_step_anon_after_release_delta = delta("session_loaded", "heatmap", "anon_after_release"),
    da_step_live_r_heap_delta = delta("session_loaded", "heatmap", "live_r_heap"),
    holders_released_anon_after_release_delta = delta("session_loaded", "holders_released", "anon_after_release"),
    holders_released_live_r_heap_delta = delta("session_loaded", "holders_released", "live_r_heap"),
    heatmap_anon_raw_minus_after_release = first("heatmap")$anon_raw - first("heatmap")$anon_after_release)
jsonlite::write_json(record, record_path,
    auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
cat(sprintf("RESULT DA step: raw +%.3f GiB, after release +%.3f GiB, live R heap +%.3f GiB; after releasing every holder: +%.3f GiB (live R heap +%.3f GiB)\n",
    gib(record$result$da_step_anon_raw_delta), gib(record$result$da_step_anon_after_release_delta),
    gib(record$result$da_step_live_r_heap_delta), gib(record$result$holders_released_anon_after_release_delta),
    gib(record$result$holders_released_live_r_heap_delta)))
