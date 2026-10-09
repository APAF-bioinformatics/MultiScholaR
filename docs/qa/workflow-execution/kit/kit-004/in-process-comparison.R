# KIT-004 in-process comparison (not a counted member).
#
# Rebuilds the `protein_replicate_filtered` state of the KIT-001 complete member (screen corpus v2,
# kit_001 control) into a fresh native store and runs normalisation + RUV, correlation filtering and
# differential abundance through the module-edge entry points with the frozen journey settings.
# Every state, product, GUI holder and publication file is fingerprinted with the volatile parts masked
# (store request ids, plot objects, timestamps, widget element ids), so two runs of the same code agree
# and a code change that alters any product, checkpoint, file or holder shows up by key.
#
#   Rscript in-process-comparison.R run <label> [out_root]
#   Rscript in-process-comparison.R compare <dir-a> <dir-b>
#
# Run from the repository root with the proof library on R_LIBS (the kit README check recipe).

args <- commandArgs(trailingOnly = TRUE)
mode <- args[[1L]]
PROOF <- "/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5"
MEMBER <- file.path(PROOF, "kit/runtime/kit-lifecycle-complete/member/project/dia_baseline_smoke")
OUT_ROOT <- if (length(args) >= 3L && mode == "run") args[[3L]] else if (length(args) >= 2L && mode == "hit") args[[2L]] else path.expand("~/.cache/multischolar-kit-004")

mask_value <- function(value, store, root, depth = 0L) {
    if (depth > 48L) return("<depth>")
    if (is.null(value)) return(NULL)
    if (is.function(value)) return("<function>")
    if (is.environment(value)) return("<environment>")
    if (inherits(value, c("gg", "ggplot", "grob", "gtable", "htmlwidget", "FilteringProgress", "pheatmap",
        "recordedplot"))) return(paste0("<", class(value)[[1L]], ">"))
    if (isS4(value)) return(list(.s4 = class(value)[[1L]], slots = mask_value(stats::setNames(lapply(
        methods::slotNames(value), \(s) methods::slot(value, s)), methods::slotNames(value)), store, root, depth + 1L)))
    if (is.list(value) && !is.data.frame(value) && isTRUE(isSavedWorkflowReference(value))) {
        if (identical(value$payload$codec, "parquet_duckdb") || is.null(store)) {
            return(list(.ref = if (is.null(store)) value$kind else "table", sha256 = value$payload$sha256,
                bytes = value$payload$bytes))
        }
        hydrated <- hydrateSavedWorkflowEntity(store, value)
        return(list(.ref = value$kind, sha256 = value$payload$sha256,
            value = mask_value(hydrated, store, root, depth + 1L)))
    }
    if (is.data.frame(value)) {
        masked <- value
        for (column in names(masked)) if (is.character(masked[[column]]))
            masked[[column]] <- gsub(root, "<root>", masked[[column]], fixed = TRUE)
        return(list(.df = digest::digest(masked, algo = "sha256"), dim = dim(value), names = names(value)))
    }
    if (is.list(value)) {
        drop <- names(value) %in% c("timestamp", "mtime")
        value <- value[!drop]
        return(lapply(value, mask_value, store = store, root = root, depth = depth + 1L))
    }
    if (is.character(value)) return(gsub(root, "<root>", value, fixed = TRUE))
    value
}

fingerprint_files <- function(project) {
    files <- list.files(project, recursive = TRUE, full.names = TRUE, all.files = TRUE)
    files <- files[!grepl("/[.]workflow-entities/|/[.]dia-[^/]*/|[.]log$", files)]
    files <- files[!dir.exists(files)]
    relative <- substring(files, nchar(project) + 2L)
    hash <- function(path) {
        lower <- tolower(path)
        if (grepl("[.]pdf$", lower)) {
            raw <- readBin(path, "raw", file.size(path))
            text <- rawToChar(raw[raw != as.raw(0L)], multiple = FALSE)
            text <- gsub("/(CreationDate|ModDate) *\\([^)]*\\)", "", text, useBytes = TRUE)
            return(c(digest::digest(text, algo = "sha256", serialize = FALSE), "pdf_dates_stripped"))
        }
        if (grepl("[.]rds$", lower)) {
            value <- tryCatch(readRDS(path), error = function(e) NULL)
            if (is.null(value)) return(c(digest::digest(file = path, algo = "sha256"), "bytes"))
            return(c(digest::digest(mask_value(value, NULL, project), algo = "sha256"), "rds_masked"))
        }
        if (grepl("[.]html$", lower)) {
            text <- readLines(path, warn = FALSE)
            text <- gsub("htmlwidget-[A-Za-z0-9]+", "htmlwidget-X", text)
            return(c(digest::digest(text, algo = "sha256", serialize = FALSE), "html_widget_ids_masked"))
        }
        if (grepl("[.]xlsx$", lower) && requireNamespace("openxlsx", quietly = TRUE)) {
            sheets <- openxlsx::getSheetNames(path)
            tables <- lapply(sheets, function(sheet) openxlsx::read.xlsx(path, sheet = sheet))
            return(c(digest::digest(tables, algo = "sha256"), "xlsx_sheets"))
        }
        c(digest::digest(file = path, algo = "sha256"), "bytes")
    }
    hashed <- vapply(files, hash, character(2))
    data.frame(path = relative, bytes = file.size(files), sha256 = hashed[1L, ], method = hashed[2L, ],
        stringsAsFactors = FALSE, row.names = NULL)
}


# Canonical form of a masked tree: store entity ids (digests of descriptors that carry the random request
# path) and rds payload hashes are dropped, hex-keyed lists are ordered by content, timestamped and
# request-id path segments are replaced.
canonical <- function(value) {
    if (is.list(value) && !is.data.frame(value)) {
        if (!is.null(value$.ref) && !identical(value$.ref, "table")) value$sha256 <- NULL
        value <- lapply(value, canonical)
        keys <- names(value)
        if (length(value) && !is.null(keys) && all(grepl("^[0-9a-f]{64}$", keys))) {
            value <- unname(value)
            value <- value[order(vapply(value, function(item) digest::digest(item, algo = "sha256"), character(1)))]
        }
        return(value)
    }
    if (is.character(value)) {
        value <- gsub("[0-9]{8}_[0-9]{6}", "<time>", value)
        value <- gsub("[0-9a-f]{32}", "<id>", value)
    }
    value
}

canonical_files <- function(files) {
    files <- files[!grepl("/[.]workflow-export/", files$path), ]
    files$path <- canonical(files$path)
    files[order(files$path), ]
}

fingerprint_store <- function(manager, holders, out) {
    store <- manager$.context$store
    snapshot <- manager$.context$snapshot
    mask <- function(value) mask_value(value, store, out)
    states <- withDIAStateScope(manager, function() stats::setNames(lapply(names(snapshot$states), function(name)
        mask(getDIAStateReference(manager, name))), names(snapshot$states)))
    products <- withDIAStateScope(manager, function() lapply(snapshot$products, function(reference) {
        value <- hydrateSavedWorkflowEntity(store, reference)
        list(product_type = value$product_type, factors = reference$metadata$factors, value = mask(value$value))
    }))
    products <- products[order(vapply(products, function(p) paste(p$product_type, digest::digest(p$value)), character(1)))]
    names(products) <- NULL
    da_results <- withDIAStateScope(manager, function() mask(snapshot$da_results))
    gui_references <- withDIAStateScope(manager, function() mask(snapshot$gui_references))
    holders_masked <- withDIAStateScope(manager, function() lapply(holders, function(owner)
        stats::setNames(lapply(names(owner), function(field) mask(owner[[field]])), names(owner))))
    ownership <- diaStateOwnership(manager)
    list(states = states, products = products, da_results = da_results, gui_references = gui_references,
        holders = holders_masked, ownership = ownership[c("full_s4_reads", "working_records")])
}

run <- function(label) {
    suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE))
    out <- file.path(OUT_ROOT, label)
    stopifnot(!dir.exists(out))
    dir.create(out, recursive = TRUE)
    timings <- list()
    started <- proc.time()[["elapsed"]]
    # 1. The source state, materialised from the immutable member store (read-only).
    source_store <- list.files(MEMBER, "^[.]dia-native-", all.files = TRUE, full.names = TRUE, recursive = TRUE,
        include.dirs = TRUE)
    source_store <- source_store[dir.exists(source_store)][[1L]]
    identity <- readSavedStoreIdentity(file.path(source_store, ".workflow-entities", "store.json"))
    head_before <- digest::digest(file = file.path(source_store, ".workflow-entities", "dia-state-head.rds"), algo = "sha256")
    source_manager <- newDIAStateAdapter(newWorkflowEntityStore(source_store, identity$project_id, identity$workflow_id,
        max_cache_bytes = 0), reopen = TRUE)
    record <- withDIAStateScope(source_manager, function() {
        reference <- getDIAStateReference(source_manager, "protein_replicate_filtered")
        materializeDIAStateRecord(source_manager$.context, reference)
    })
    closeDIAStateAdapter(source_manager)
    stopifnot(identical(head_before, digest::digest(file = file.path(source_store, ".workflow-entities",
        "dia-state-head.rds"), algo = "sha256")))
    object <- record$data
    stopifnot(methods::is(object, "ProteinQuantitativeData"))
    timings$source <- proc.time()[["elapsed"]] - started
    # 2. A fresh project and native store holding only that state.
    project <- file.path(out, "project")
    dir.create(project)
    directories <- setupDirectories(project, "proteomics")
    paths <- directories$proteomics
    directories$proteomics_dia_baseline_smoke <- paths
    state_root <- file.path(paths$source_dir, ".dia-native-kit004")
    dir.create(state_root)
    manager <- newDIAStateAdapter(newWorkflowEntityStore(state_root, "kit-004", "DIA"), WorkflowState$new(audit_enabled = FALSE))
    manager$.context$reference_states <- TRUE
    manager$.context$experiment_paths <- paths
    manager$.context$reclaim_enabled <- TRUE
    manager$saveState("protein_replicate_filtered", object, record$config, record$description)
    scripts <- file.path(MEMBER, "scripts", "proteomics")
    sidecars <- c("aa_seq_tbl_final", "uniprot_dat_cln", "fasta_metadata")
    annotation <- stats::setNames(lapply(sidecars, function(name) readRDS(file.path(scripts, paste0(name, ".RDS")))), sidecars)
    contrasts <- readProtDesignImportedContrasts(file.path(scripts, "contrast_strings.tab"))
    for (name in c("config.ini", "contrast_strings.tab", "design_matrix.tab", paste0(sidecars, ".RDS")))
        stopifnot(file.copy(file.path(scripts, name), paths$source_dir, overwrite = TRUE))
    config <- object@args
    object <- NULL
    withr::local_options(multischolar.test_mode = FALSE, multischolar.dia_native_sessions = TRUE,
        multischolar.capture_test_checkpoints = FALSE)
    withr::local_envvar(MULTISCHOLAR_TEST_MODE = "false")
    rlang::local_bindings(config_list = config, project_dirs = directories,
        filtering_progress = methods::new("FilteringProgress"), .env = .GlobalEnv)
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    settings <- list(norm_method = "cyclicloess", ruv_mode = "automatic", auto_percentage_min = 1,
        auto_percentage_max = 20, separation_metric = "max_difference", k_penalty_weight = 0.5,
        adaptive_k_penalty = TRUE, color_variable = "group", shape_variable = NULL,
        min_pearson_correlation_threshold = 0.5)
    aesthetics <- function() list(color_var = "group", shape_var = NULL)
    holders <- list()
    fingerprints <- NULL
    shiny::testServer(function(input, output, session) {
        design <- readDIAExecutionComponents(manager, fields = "design_matrix")@design_matrix
        workflow <- shiny::reactiveValues(state_manager = manager, design_matrix = design, config_list = config,
            contrasts_tbl = contrasts, data_type = "peptide", data_format = "diann",
            aa_seq_tbl_final = annotation$aa_seq_tbl_final, uniprot_dat_cln = annotation$uniprot_dat_cln,
            fasta_metadata = annotation$fasta_metadata)
        attachDIAStateSession(manager, session, workflow)
        norm <- shiny::reactiveValues()
        da <- shiny::reactiveValues()
        da_server_run_analysis_handler(input, output, session, session$ns, da, workflow, paths)
        session$userData$phases <- list(
            normalisation = function() shiny::isolate({
                generateProtNormPreNormalizationQcArtifacts(manager, paths$protein_qc_dir, aesthetics())
                post <- function(reference) {
                    norm$qc_plot_paths <- generateProtNormPostNormalizationQcArtifacts(reference,
                        norm$publication_paths$protein_qc_dir %||% paths$protein_qc_dir, aesthetics(), norm$qc_plot_paths)
                }
                ruv <- function(reference) {
                    norm$qc_plot_paths <- generateProtNormRuvCorrectedQcArtifacts(reference,
                        norm$publication_paths$protein_qc_dir %||% paths$protein_qc_dir, aesthetics(), norm$qc_plot_paths)
                }
                stopifnot(isTRUE(runProtNormNormalizationWorkflow(settings, workflow, norm, paths, "proteomics",
                    "dia_baseline_smoke", checkMemoryUsageFn = function(...) NULL, generatePostNormalizationQcFn = post,
                    generateRuvCorrectedQcFn = ruv, getRuvGroupingVariableFn = function() "group",
                    showNotificationFn = function(...) NULL, handleErrorFn = function(error, ...) stop(error))))
            }),
            correlation = function() shiny::isolate({
                final <- runProtNormApplyCorrelationObserver(settings, output, workflow, norm, paths, "proteomics",
                    "dia_baseline_smoke", function() "group", aesthetics, checkMemoryUsageFn = function(...) NULL,
                    handleCorrelationErrorFn = function(error, ...) stop(error))
                stopifnot(isSavedWorkflowReference(final))
            }),
            da_prepare = function() shiny::isolate({
                da$current_s4_object <- getDIAStateReference(manager)
                da$contrasts_available <- contrasts$contrasts
            }),
            collect = function() shiny::isolate({
                holders <<- list(norm = shiny::reactiveValuesToList(norm), da = shiny::reactiveValuesToList(da),
                    workflow = shiny::reactiveValuesToList(workflow)[c("ruv_optimization_result", "protein_counts",
                        "da_ui_params", "tab_status", "da_analysis_results_list")])
            }))
    }, {
        for (phase in c("normalisation", "correlation", "da_prepare")) {
            t0 <- proc.time()[["elapsed"]]
            session$userData$phases[[phase]]()
            session$flushReact()
            timings[[phase]] <<- proc.time()[["elapsed"]] - t0
        }
        t0 <- proc.time()[["elapsed"]]
        session$setInputs(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0, run_da_analysis = 1L)
        session$flushReact()
        timings$da <<- proc.time()[["elapsed"]] - t0
        stopifnot(isTRUE(shiny::isolate(da$analysis_complete)))
        session$userData$phases$collect()
        fingerprints <<- fingerprint_store(manager, holders, out)
    })
    files <- canonical_files(fingerprint_files(project))
    fingerprints <- canonical(fingerprints)
    states <- fingerprints$states
    products <- fingerprints$products
    da_results <- fingerprints$da_results
    gui_references <- fingerprints$gui_references
    holders_masked <- fingerprints$holders
    result <- list(label = label, head_commit = system("git rev-parse --short HEAD", intern = TRUE),
        source_member_head_sha256 = head_before, source_state = "protein_replicate_filtered",
        settings = settings, da = list(formula = "~ 0 + group", qvalue = 0.05, lfc = 0),
        states = states, products = products, da_results = da_results, gui_references = gui_references,
        holders = holders_masked, files = files, timings = timings,
        ownership = fingerprints$ownership, elapsed = proc.time()[["elapsed"]] - started)
    saveRDS(result, file.path(out, "fingerprints.rds"))
    digests <- list(label = label, head_commit = result$head_commit, source_member_head_sha256 = head_before,
        states = lapply(states, digest::digest, algo = "sha256"),
        products = lapply(products, function(p) list(type = p$product_type, sha256 = digest::digest(p$value, algo = "sha256"))),
        da_results = digest::digest(da_results, algo = "sha256"),
        gui_references = lapply(gui_references, digest::digest, algo = "sha256"),
        holders = lapply(holders_masked, function(owner) lapply(owner, digest::digest, algo = "sha256")),
        files = stats::setNames(as.list(files$sha256), files$path), timings = timings, elapsed = result$elapsed)
    jsonlite::write_json(digests, file.path(out, "digests.json"), auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA)
    cat(sprintf("%s: %d states, %d products, %d files, %.0f s\n", label, length(states), length(products), nrow(files), result$elapsed))
}

diff_paths <- function(a, b, path = "") {
    if (identical(a, b)) return(character())
    if (is.list(a) && is.list(b) && !is.data.frame(a) && !is.data.frame(b)) {
        keys <- union(names(a) %||% as.character(seq_along(a)), names(b) %||% as.character(seq_along(b)))
        if (is.null(names(a)) || is.null(names(b))) {
            if (length(a) != length(b)) return(paste0(path, " (length ", length(a), " vs ", length(b), ")"))
            return(unlist(lapply(seq_along(a), function(i) diff_paths(a[[i]], b[[i]], paste0(path, "[[", i, "]]")))))
        }
        return(unlist(lapply(keys, function(key) {
            if (!key %in% names(a) || !key %in% names(b)) return(paste0(path, "$", key, " (missing on one side)"))
            diff_paths(a[[key]], b[[key]], paste0(path, "$", key))
        })))
    }
    path
}

compare <- function(dir_a, dir_b) {
    a <- canonical(readRDS(file.path(dir_a, "fingerprints.rds")))
    b <- canonical(readRDS(file.path(dir_b, "fingerprints.rds")))
    a$files <- canonical_files(a$files)
    b$files <- canonical_files(b$files)
    differences <- character()
    for (section in c("states", "products", "da_results", "gui_references", "holders")) {
        found <- diff_paths(a[[section]], b[[section]], section)
        differences <- c(differences, found)
    }
    files <- merge(a$files[c("path", "sha256")], b$files[c("path", "sha256")], by = "path", all = TRUE, suffixes = c("_a", "_b"))
    changed <- files[is.na(files$sha256_a) | is.na(files$sha256_b) | files$sha256_a != files$sha256_b, ]
    differences <- c(differences, if (nrow(changed)) paste0("files$", changed$path))
    cat(sprintf("compare %s (%s) vs %s (%s): %d states, %d products, %d files\n", a$label, a$head_commit, b$label,
        b$head_commit, length(a$states), length(a$products), nrow(a$files)))
    if (length(differences)) {
        cat("DIFFERENCES:\n"); cat(paste0("  ", differences, "\n"), sep = "")
        quit(status = 1L)
    }
    cat("IDENTICAL (states, products, DA manifest, GUI references, holders, files)\n")
}


# Identity continuity (audit B4): every product the KIT-001 complete member saved on the pre-ticket
# build must hit under the runner, with the member's own recorded parameters, without any kernel call.
#   Rscript in-process-comparison.R hit [out_root]
hit <- function() {
    suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE))
    out <- file.path(OUT_ROOT, "hit")
    stopifnot(!dir.exists(out))
    dir.create(out, recursive = TRUE)
    source_store <- list.files(MEMBER, "^[.]dia-native-", all.files = TRUE, full.names = TRUE, recursive = TRUE,
        include.dirs = TRUE)
    source_store <- file.path(source_store[dir.exists(source_store)][[1L]], ".workflow-entities")
    # Copy of the member store (payloads hard-linked, never written): the member stays immutable.
    target <- file.path(out, ".workflow-entities")
    dir.create(target)
    for (name in c("store.json", "dia-state-head.rds")) stopifnot(file.copy(file.path(source_store, name), target))
    files <- list.files(file.path(source_store, "entities"), recursive = TRUE, full.names = TRUE)
    for (path in files[!dir.exists(files)]) {
        to <- file.path(target, substring(path, nchar(source_store) + 2L))
        dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
        if (grepl("/payload[.](parquet|rds|qs2)$", path)) stopifnot(file.link(path, to)) else stopifnot(file.copy(path, to))
    }
    identity <- readSavedStoreIdentity(file.path(target, "store.json"))
    manager <- newDIAStateAdapter(newWorkflowEntityStore(out, identity$project_id, identity$workflow_id), reopen = TRUE)
    on.exit(closeDIAStateAdapter(manager), add = TRUE)
    manager$.context$reference_states <- TRUE
    kernel_called <- function(...) stop("a kernel ran on an identity hit")
    report <- testthat::with_mocked_bindings(withDIAStateScope(manager, function() {
        source <- getDIAStateReference(manager, "protein_replicate_filtered")
        ruv_state <- readDIAStateMetadata(manager, getDIAStateReference(manager, "ruv_corrected"))$data
        expected_product <- ruv_state$products$normalization_science
        # The software identity is the installed package binary's digest; the member ran on the installed
        # kit_001 build and this script runs from the working tree, so the member's software identity is
        # substituted. What is tested is the identity scheme (operation names, parameters, slots, selectors).
        manager$.context$software <- expected_product$identity$software
        saved <- hydrateSavedWorkflowEntity(manager$.context$store, expected_product)$value
        normalization <- runSavedOperation(manager, diaNormalizationOperation(), source, saved$parameters)
        stopifnot(identical(normalization$status, "hit"), identical(normalization$reference, expected_product),
            identical(normalization$value$final, saved$final))
        final <- saved$final
        pairs <- runSavedOperation(manager, diaPearsonPairsOperation(), final, list(group = "group", tech_rep_remove_regex = "pool"))
        stopifnot(identical(pairs$status, "hit"))
        selection <- runSavedOperation(manager, diaPearsonSelectionOperation(), final,
            list(pairs = pairs$reference$entity_id, threshold = 0.5, correlations = pairs$value), dependencies = list(pairs$reference))
        correlation_state <- readDIAStateMetadata(manager, getDIAStateReference(manager, "correlation_filtered"))$data
        stopifnot(identical(selection$status, "hit"),
            identical(readDIAStateMetadata(manager, selection$value)$data$components, correlation_state$components))
        manifest <- hydrateSavedWorkflowEntity(manager$.context$store, manager$.context$snapshot$da_results)
        spec <- readDIAStateMetadata(manager, manifest$source)$data
        contrasts <- manifest$contrasts
        da_hits <- vapply(seq_len(nrow(contrasts)), function(i) {
            run <- runSavedOperation(manager, diaContrastOperation(), manifest$source, c(manifest$parameters,
                list(contrast = contrasts$full_format[[i]], annotations = spec$components$protein_id_table$entity_id,
                    name = contrasts$contrasts[[i]], comparison = contrasts$friendly_names[[i]])))
            identical(run$status, "hit") && identical(run$value, manifest$results[[contrasts$contrasts[[i]]]])
        }, logical(1))
        stopifnot(all(da_hits))
        list(normalization = expected_product$entity_id, pairs = pairs$reference$entity_id,
            selection = selection$reference$entity_id, da_contrasts = nrow(contrasts), da_hits = sum(da_hits),
            parameters = saved$parameters, software_identity = "substituted from the member's product (installed kit_001 build)",
            full_s4_reads = manager$.context$full_s4_reads)
    }), normaliseBetweenSamples = kernel_called, pearsonCorForSamplePairs = kernel_called,
        filterSamplesByProteinCorrelationThreshold = kernel_called, differentialAbundanceAnalysis = kernel_called)
    report$head_commit <- system("git rev-parse --short HEAD", intern = TRUE)
    report$member_head_sha256 <- digest::digest(file = file.path(source_store, "dia-state-head.rds"), algo = "sha256")
    jsonlite::write_json(report, file.path(out, "identity-hits.json"), auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA)
    cat(sprintf("identity hits: normalization, pearson pairs, pearson selection, %d/%d DA contrasts; full S4 reads %d\n",
        report$da_hits, report$da_contrasts, report$full_s4_reads))
}

if (mode == "run") run(args[[2L]]) else if (mode == "compare") compare(args[[2L]], args[[3L]]) else if (mode == "hit") hit() else stop("unknown mode")
