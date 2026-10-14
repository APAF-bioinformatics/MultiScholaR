# Harness-only instrumentation. Never sourced by the production package.

#' Record one durable, payload-free diagnostic boundary without forcing GC.
#' @param context Diagnostic owner containing paths and scalar counters only.
#' @param stage Named function or lifecycle boundary.
#' @param detail Bounded metadata, never assay values or environments.
#' @return NULL invisibly.
#' @noRd
checkpointGuiDiagnosticEvent <- function(context, stage, detail = list()) {
    if (isTRUE(context$recording_limited)) return(invisible(NULL))
    path <- file.path(context$directory, "application-phases.jsonl")
    if (file.exists(path) && file.size(path) >= context$max_log_bytes) {
        context$recording_limited <- TRUE
        stage <- "diagnostic:recording_limit"
        detail <- list(limit_bytes = context$max_log_bytes, coverage_complete = FALSE)
    }
    status <- readLines("/proc/self/status", warn = FALSE)
    fields <- c("VmRSS", "VmHWM", "RssAnon", "RssFile", "RssShmem")
    memory <- stats::setNames(lapply(fields, function(field) {
        line <- grep(paste0("^", field, ":"), status, value = TRUE)
        if (length(line) != 1L) return(NA_real_)
        as.numeric(strsplit(trimws(sub("^[^:]+:", "", line)), " +")[[1L]][1L]) * 1024
    }), fields)
    arrow_pool <- NULL
    if ("arrow" %in% loadedNamespaces()) {
        pool <- arrow::default_memory_pool()
        arrow_pool <- list(backend = pool$backend_name,
            allocated_bytes = pool$bytes_allocated, peak_bytes = pool$max_memory)
    }
    record <- list(time = as.numeric(Sys.time()), pid = Sys.getpid(), stage = stage,
        cpu = unname(proc.time()[1:2]), memory = memory, arrow_pool = arrow_pool,
        sql_queries = context$sql_queries, detail = detail)
    connection <- file(path, "a")
    on.exit(close(connection), add = TRUE)
    writeLines(jsonlite::toJSON(record, auto_unbox = TRUE, null = "null", digits = NA), connection)
    flush(connection)
    flush(context$console)
    invisible(NULL)
}

#' Enable bounded glibc attribution in an explicitly selected diagnostic process.
#' @param context Existing diagnostic context; never a production application owner.
#' @param arena_details Include bounded all-arena and Linux mapping summaries.
#' @return NULL invisibly.
#' @noRd
enableCheckpointAllocatorProbe <- function(context, arena_details = FALSE) {
    native <- new.env(parent = baseenv())
    Rcpp::cppFunction(code = '
Rcpp::NumericVector allocatorSnapshot(bool trim) {
    int released = trim ? malloc_trim(0) : -1;
    struct mallinfo2 info = mallinfo2();
    return Rcpp::NumericVector::create(
        Rcpp::Named("arena_bytes") = (double) info.arena,
        Rcpp::Named("in_use_arena_bytes") = (double) info.uordblks,
        Rcpp::Named("free_arena_bytes") = (double) info.fordblks,
        Rcpp::Named("mmap_bytes") = (double) info.hblkhd,
        Rcpp::Named("main_top_free_bytes") = (double) info.keepcost,
        Rcpp::Named("free_chunk_count") = (double) info.ordblks,
        Rcpp::Named("trim_return") = released);
}', includes = "#include <malloc.h>", env = native)
    context$allocator_snapshot <- native$allocatorSnapshot
    if (isTRUE(arena_details)) {
        Rcpp::cppFunction(code = '
std::string allocatorXml() {
    std::vector<char> buffer(1048576, 0);
    FILE* stream = fmemopen(buffer.data(), buffer.size(), "w");
    if (!stream) Rcpp::stop("cannot open bounded allocator buffer");
    int status = malloc_info(0, stream);
    int flushed = fflush(stream);
    long count = ftell(stream);
    int closed = fclose(stream);
    if (status || flushed || closed || count < 0 || count >= 1048576)
        Rcpp::stop("allocator metadata failed or exceeded 1 MiB");
    return std::string(buffer.data(), (size_t) count);
}', includes = "#include <malloc.h>\n#include <stdio.h>", env = native)
        context$allocator_xml <- native$allocatorXml
    }
    context$allocator_pending <- FALSE
    context$allocator_started <- FALSE
    context$allocator_trimmed <- FALSE
    checkpointGuiDiagnosticEvent(context, "allocator:enabled", list(timing_eligible = FALSE))
    invisible(NULL)
}

#' Summarize glibc arenas without retaining its free-block XML inventory.
#' @param xml Bounded malloc_info XML from the diagnostic native counter.
#' @return Per-arena virtual allocation/free bytes and aggregate cross-check.
#' @noRd
checkpointAllocatorArenas <- function(xml) {
    stopifnot(nchar(xml, type = "bytes") <= 1048576)
    doc <- xml2::read_xml(xml)
    heaps <- xml2::xml_find_all(doc, "./heap")
    stopifnot(length(heaps) > 0L, length(heaps) <= 1024L)
    bytes <- function(node, path) as.numeric(xml2::xml_attr(xml2::xml_find_all(node, path), "size"))
    arenas <- lapply(heaps, function(heap) list(
        id = as.integer(xml2::xml_attr(heap, "nr")),
        system_bytes = bytes(heap, "./system[@type='current']"),
        free_bytes = sum(bytes(heap, "./total[@type='rest' or @type='fast']"))))
    system <- bytes(doc, "./system[@type='current']")
    stopifnot(sum(vapply(arenas, `[[`, numeric(1), "system_bytes")) == system)
    list(count = length(arenas), system_bytes = system,
        free_bytes = sum(bytes(doc, "./total[@type='rest' or @type='fast']")), arenas = arenas)
}

#' Aggregate Linux memory mappings without recording addresses or file names.
#' @return RSS/anonymous/private-dirty byte totals in four mapping categories.
#' @noRd
checkpointAllocatorMappings <- function() {
    categories <- c("main_heap", "anonymous_mapping", "stack", "other_mapping")
    fields <- c(Rss = "rss_bytes", Anonymous = "anonymous_bytes", Private_Dirty = "private_dirty_bytes")
    totals <- matrix(0, length(categories), length(fields), dimnames = list(categories, unname(fields)))
    input <- file("/proc/self/smaps", "r")
    on.exit(close(input), add = TRUE)
    category <- "other_mapping"
    seen <- 0L
    repeat {
        lines <- readLines(input, n = 256L, warn = FALSE)
        if (!length(lines)) break
        seen <- seen + length(lines)
        stopifnot(seen <= 100000L)
        for (line in lines) {
            if (grepl("^[[:xdigit:]]+-[[:xdigit:]]+ ", line)) {
                parts <- strsplit(trimws(line), " +")[[1L]]
                category <- if (grepl("[heap]", line, fixed = TRUE)) "main_heap" else
                    if (grepl("[stack", line, fixed = TRUE)) "stack" else
                        if (length(parts) <= 5L || grepl("[anon", line, fixed = TRUE))
                            "anonymous_mapping" else "other_mapping"
            } else {
                name <- sub(":.*", "", line)
                if (!name %in% names(fields)) next
                value <- as.numeric(strsplit(trimws(line), " +")[[1L]][[2L]]) * 1024
                totals[category, fields[[name]]] <- totals[category, fields[[name]]] + value
            }
        }
    }
    stats::setNames(lapply(categories, function(name) as.list(totals[name, ])), categories)
}

#' Measure GC and allocator residency without retaining the observed call frame.
#' @param context Explicit allocator diagnostic context.
#' @param label Completed-action or figure boundary.
#' @param collect Whether to measure live R heap after full GC.
#' @param trim Whether to release glibc free pages in this diagnostic process only.
#' @return NULL invisibly.
#' @noRd
checkpointAllocatorBoundary <- function(context, label, collect = FALSE, trim = FALSE) {
    if (!is.function(context$allocator_snapshot)) return(invisible(NULL))
    live <- NULL
    if (collect) {
        cells <- gc(full = TRUE)
        live <- sum(cells[c("Ncells", "Vcells"), "used"] * c(56, 8))
    }
    native <- as.list(context$allocator_snapshot(trim))
    arenas <- maps <- NULL
    if (is.function(context$allocator_xml)) {
        arenas <- checkpointAllocatorArenas(context$allocator_xml())
        maps <- checkpointAllocatorMappings()
    }
    checkpointGuiDiagnosticEvent(context, "allocator:boundary",
        list(label = label, live_R_bytes = live, collected = collect, trimmed = trim,
            glibc = native, arenas = arenas, mappings = maps))
    invisible(NULL)
}

#' Describe named owners without evaluating reactive getters or manager bindings.
#' @param values Already available ordinary list of roots.
#' @return Shallow class/address/shape metadata; bytes are deliberately unmeasured.
#' @noRd
checkpointGuiOwnerSizes <- function(values) {
    lapply(values, function(value) {
        list(class = class(value), address = lobstr::obj_addr(value),
            bytes = NA_real_,
            shape = if (is.data.frame(value) || is.matrix(value)) dim(value) else NULL)
    })
}

#' Observe completed-action live heap separately from resident process memory.
#' @param context Diagnostic context; collection is opt-in and never a timed run.
#' @param label Completed action/state label.
#' @return Post-collection live bytes, not RSS or cumulative allocations.
#' @noRd
checkpointGuiCompletedHeap <- function(context, label) {
    if (!isTRUE(context$collect_heap)) return(invisible(NULL))
    checkpointGuiDiagnosticEvent(context, "heap:before_collection", list(label = label))
    cells <- gc(full = TRUE)
    live <- sum(cells[c("Ncells", "Vcells"), "used"] * c(56, 8))
    checkpointGuiDiagnosticEvent(context, "heap:after_collection",
        list(label = label, live_heap_bytes = live, gc_cells = unname(cells[, "used"]),
            note = "64-bit R cell accounting after full collection; not a timed member"))
    invisible(live)
}

#' Inspect only named ordinary owners while their producer frame is already live.
#' @param frame Current helper frame; never retained after observation.
#' @return Bounded shallow metadata, not scientific objects.
#' @noRd
checkpointGuiFrameOwners <- function(frame) {
    fields <- intersect(c("normData", "da_data", "workflowData", "sessionData"),
        ls(frame, all.names = TRUE))
    stats::setNames(lapply(fields, function(field) {
        value <- get(field, frame, inherits = FALSE)
        if (inherits(value, "reactivevalues")) {
            value <- shiny::isolate(shiny::reactiveValuesToList(value, all.names = TRUE))
        } else if (is.environment(value)) {
            names <- ls(value, all.names = TRUE)
            names <- names[!vapply(names, bindingIsActive, logical(1), env = value)]
            value <- mget(utils::head(names, 64L), value, inherits = FALSE)
        }
        if (is.list(value)) checkpointGuiOwnerSizes(utils::head(value, 64L)) else
            checkpointGuiOwnerSizes(list(value = value))
    }), fields)
}

#' Attach post-delivery owner observations to the existing session only.
#' @param context Diagnostic context; never stores workflow payloads.
#' @param frame Completed initializer call frame.
#' @return NULL invisibly.
#' @noRd
checkpointGuiObserveOwners <- function(context, frame) {
    session <- get("session", frame, inherits = FALSE)
    workflow <- get("workflow_data", frame, inherits = FALSE)
    session$onFlushed(function() {
        values <- shiny::isolate(shiny::reactiveValuesToList(workflow, all.names = TRUE))
        manager <- session$userData$dia_capacity_manager
        detail <- list(owners = checkpointGuiOwnerSizes(values),
            note = "shallow metadata only; no object.size, forced GC or hydration; not a heap census")
        if (inherits(manager, "DIAStateAdapter")) {
            owner <- manager$.context
            detail$state <- owner$snapshot$current_state
            detail$working_ids <- ls(owner$working, all.names = TRUE)
            detail$borrowed_worker_groups <- names(owner$borrowed_workers)
            detail$owned_worker_groups <- names(owner$owned_workers)
        }
        checkpointGuiDiagnosticEvent(context, "flush:owners", detail)
        if (isTRUE(context$design_memory)) {
            label <- shiny::isolate(session$rootScope()$input$.dia_memory_probe)
            if (length(label) == 1L && label %in% c("import", "design_setup") &&
                !label %in% context$design_memory_done &&
                inherits(manager, "DIAStateAdapter") && manager$.context$scope_depth == 0L &&
                manager$.context$store$open_resources == 0L &&
                is.null(manager$.context$store$table_action)) {
                values <- NULL
                checkpointAllocatorBoundary(context, paste0(label, ":observed"))
                if (label == "design_setup") {
                    checkpointAllocatorBoundary(context, paste0(label, ":gc"), collect = TRUE)
                    checkpointAllocatorBoundary(context, paste0(label, ":trim"), trim = TRUE)
                }
                context$design_memory_done <- c(context$design_memory_done, label)
                jsonlite::write_json(list(status = if (isTRUE(context$recording_limited)) "incomplete" else "complete", label = label),
                    file.path(context$directory, paste0("memory-boundary-", label, ".json")),
                    auto_unbox = TRUE)
            }
        }
        if (isTRUE(context$allocator_pending)) {
            context$allocator_pending <- FALSE
            values <- NULL
            checkpointAllocatorBoundary(context, paste0("delivered:", detail$state))
            checkpointAllocatorBoundary(context, paste0("delivered_gc:", detail$state), collect = TRUE)
            if (identical(detail$state, "qvalue_filtered") && !isTRUE(context$allocator_trimmed)) {
                checkpointAllocatorBoundary(context, "qvalue_delivered_trim", collect = TRUE, trim = TRUE)
                context$allocator_trimmed <- TRUE
            }
        }
        if (length(context$heap_pending)) {
            label <- paste(context$heap_pending, collapse = ";")
            context$heap_pending <- character()
            # Do not keep the shallow observation's payload list during collection.
            values <- NULL
            checkpointGuiCompletedHeap(context, label)
        }
    }, once = FALSE)
    invisible(NULL)
}

#' Trace the installed functions used by the approved journey without replacing science.
#' @param context Diagnostic context.
#' @param namespace Installed package namespace, or isolated fixture environment.
#' @param functions Exact function names; NULL selects the approved journey.
#' @return NULL invisibly.
#' @noRd
installCheckpointGuiTraces <- function(context, namespace, functions = NULL) {
    if (is.null(functions)) functions <- c(
        "importProtImportDataByFormat", "applyProtImportResultToWorkflow",
        "processProtImportFastaData", "buildProtDesignStateCheckpoint",
        "completeProtDesignPostCheckpoint",
        "saveDIAStateValue", "resolveDIAStateValue", "saveDIAStateRecord",
        "getDIAExecutionState", "startDIAWorkerCluster", "materializeDIAStates",
        "runPeptideQvalueApplyStep", "runProteinPeptideApplyStep",
        "runPeptideRollupApplyStep", "runPeptideIntensityApplyStep",
        "runPeptideSampleApplyStep", "runPeptideReplicateApplyStep",
        "runPeptideImputationStep", "runProtNormNormalizationWorkflow",
        "buildProtNormQcPlotProduct", "buildProtNormCorrelationPlot",
        "runProtNormBetweenSamplesStep", "runProtNormExportSessionWorkflow",
        "collectProtNormExportSessionData", "saveProtNormExportArtifacts",
        "runDIAOwnedSessionExport", "hydrateSavedWorkflowEntity",
        "validateSavedWorkflowReference", "saveDIAOwnedValue",
        "writeSavedWorkflowQsPayload", "stageWorkflowLegacyExport",
        "workflowExportBoundary", "publishWorkflowLegacyExport", "resolveDIAOwnedExport",
        "restoreDIAStateSession", "runTestsContrasts",
        ".emitPeptideQcAuditRecord", ".peptideQcDigest", ".peptideQcCanonicalise",
        ".peptideQcAuditSummary", ".peptideQcStableIdentities", ".peptideQcIdentityKey",
        ".peptideQcRemovalLedger", "runProteinIqRollupApplyStep",
        "runProteinIntensityFilterApplyStep", "runProteinReplicateFilterApplyStep",
        "generateProtNormCompositeQcFigure", "generateProtNormCompositeFromFiles",
        "loadProtNormImageAsPlot", "savePlot", "runProtSummaryPublicationCopy")
    absent <- functions[!vapply(functions, exists, logical(1), envir = namespace, inherits = FALSE)]
    if (length(absent)) stop(paste("missing diagnostic functions:", paste(absent, collapse = ", ")))
    observer <- function(phase, name, frame) {
        if (name == ".peptideQcCanonicalise" &&
            !is.data.frame(get("value", frame, inherits = FALSE))) return(invisible(NULL))
        detail <- list()
        if (name %in% c("saveDIAStateValue", "saveDIAOwnedValue") && phase == "enter") {
            detail$label <- get("label", frame, inherits = FALSE)
        }
        if (name %in% c("resolveDIAStateValue", "hydrateSavedWorkflowEntity",
            "validateSavedWorkflowReference")) {
            reference <- get("reference", frame, inherits = FALSE)
            detail$entity_id <- reference$entity_id
            detail$label <- reference$identity$parameters
            detail$codec <- reference$payload$codec
            detail$payload_bytes <- reference$payload$bytes
        }
        if (name %in% c("resolveDIAStateValue", "materializeDIAStates")) {
            owner <- get("context", frame, inherits = FALSE)
            detail$export_shared_tables <- owner$export_shared_tables
            detail$scope_depth <- owner$scope_depth
            detail$working_records <- length(ls(owner$working, all.names = TRUE))
        }
        if (name == "resolveDIAStateValue") {
            detail$working_hit <- exists(reference$entity_id, owner$working, inherits = FALSE)
        }
        if (name == "materializeDIAStates") detail$share_immutable <-
            get("share_immutable", frame, inherits = FALSE)
        if (name == "workflowExportBoundary") detail$boundary <-
            get("phase", frame, inherits = FALSE)
        if (name == "writeSavedWorkflowQsPayload") detail$file_bytes <-
            unname(file.size(get("path", frame, inherits = FALSE)))
        if (name == "buildProtNormQcPlotProduct") detail$plot_type <- get("plotType", frame)
        if (name == "savePlot") detail$plot_name <- get("plot_name", frame)
        if (isTRUE(context$collect_heap) && phase == "exit" &&
            name %in% c("applyProtImportResultToWorkflow", "completeProtDesignPostCheckpoint",
                "runPeptideQvalueApplyStep", "runProteinPeptideApplyStep", "runPeptideRollupApplyStep",
                "runPeptideIntensityApplyStep", "runPeptideSampleApplyStep", "runPeptideReplicateApplyStep",
                "runPeptideImputationStep", "runProteinIqRollupApplyStep", "runProteinReplicateFilterApplyStep",
                "runProtNormNormalizationWorkflow", "runProtNormExportSessionWorkflow",
                "restoreDIAStateSession", "runTestsContrasts", "runProtSummaryPublicationCopy")) {
            context$heap_pending <- unique(c(context$heap_pending, name))
            detail$frame_owners <- checkpointGuiFrameOwners(frame)
        }
        checkpointGuiDiagnosticEvent(context, paste(phase, name, sep = ":"), detail)
        if (is.function(context$allocator_snapshot)) {
            if (phase == "enter" && name == "workflowExportBoundary" &&
                identical(detail$boundary, "before_qvalue_figures") && !context$allocator_started) {
                context$allocator_started <- TRUE
                checkpointAllocatorBoundary(context, "before_first_qvalue_figures")
                checkpointAllocatorBoundary(context, "before_first_qvalue_figures_gc", collect = TRUE)
            }
            if (phase == "enter" && name == "workflowExportBoundary" &&
                identical(detail$boundary, "before_qvalue_acknowledgement")) {
                checkpointAllocatorBoundary(context, "figures_published")
                context$allocator_pending <- TRUE
            }
            if (phase == "exit" && name %in% c("runProtNormNormalizationWorkflow",
                "runTestsContrasts", "runProtSummaryPublicationCopy")) {
                checkpointAllocatorBoundary(context, paste0("action_exit:", name))
                context$allocator_pending <- TRUE
            }
        }
        # KIT-011 removed the capacity session; the owner observation now follows the import adoption (KIT-034, Q-20).
        if (phase == "exit" && name == "applyProtImportResultToWorkflow") {
            checkpointGuiObserveOwners(context, frame)
        }
    }
    options(multischolar.gui_diagnostic_observer = observer)
    for (name in functions) {
        enter <- substitute(getOption("multischolar.gui_diagnostic_observer")("enter", NAME, environment()),
            list(NAME = name))
        leave <- substitute(getOption("multischolar.gui_diagnostic_observer")("exit", NAME, environment()),
            list(NAME = name))
        base::trace(name, tracer = enter, exit = leave, print = FALSE, where = namespace)
        context$traced <- c(context$traced, name)
    }
    context$namespace <- namespace
    if ("duckdb" %in% loadedNamespaces() && !isTRUE(context$sql_traced)) {
        context$previous_sql_observer <- getOption("multischolar.gui_sql_observer")
        options(multischolar.gui_sql_observer = function() {
            context$sql_queries <- context$sql_queries + 1
            invisible(NULL)
        })
        base::trace("dbSendQuery", signature = c("duckdb_connection", "character"),
            where = asNamespace("duckdb"), print = FALSE,
            tracer = quote(getOption("multischolar.gui_sql_observer")()))
        context$sql_traced <- TRUE
    }
    checkpointGuiDiagnosticEvent(context, "traces:installed", list(functions = functions))
    invisible(NULL)
}

#' Describe only already-bound matrices and scalar tuning values in one call.
#' @param frame Current scientific frame; it is never saved in the diagnostic owner.
#' @return Small shape/address metadata, not assays, plots or heap-size traversals.
#' @noRd
checkpointRuvPearsonShapes <- function(frame) {
    fields <- intersect(c("Y", "M", "x", "y", "mat", "cor_mat", "cln_mat", "cln_mat_4",
        "normalised_frozen_protein_matrix_filt", "theObject", "normalizedS4",
        "normalised_protein_matrix_obj", "ruvS4", "ruvCorrectedS4"), ls(frame, all.names = TRUE))
    detail <- list()
    for (field in fields) {
        value <- get(field, frame, inherits = FALSE)
        if (isS4(value) && "protein_quant_table" %in% methods::slotNames(value)) {
            value <- methods::slot(value, "protein_quant_table")
        }
        if (is.matrix(value) || is.data.frame(value)) {
            detail[[field]] <- list(class = class(value), shape = dim(value),
                address = rlang::obj_address(value),
                numeric_payload_bytes = if (is.matrix(value) && is.double(value)) length(value) * 8 else NULL)
        }
    }
    for (field in intersect(c("k", "ruv_k", "ruv_number_k", "percentage_as_neg_ctrl"),
        ls(frame, all.names = TRUE))) {
        value <- get(field, frame, inherits = FALSE)
        if (is.numeric(value) && length(value) == 1L) detail[[field]] <- value
    }
    for (field in intersect(c("ctrl", "ctl", "potentialControls", "control_genes_index"),
        ls(frame, all.names = TRUE))) {
        value <- get(field, frame, inherits = FALSE)
        detail[[field]] <- list(length = length(value),
            selected = if (is.logical(value)) sum(value, na.rm = TRUE) else NULL)
    }
    detail
}

#' Trace real RUV/Pearson methods and the exact imported numerical calls they use.
#' @param context Existing diagnostic owner.
#' @param namespace Installed MultiScholaR namespace.
#' @return NULL invisibly; all bindings/options are restored by the existing stop helper.
#' @noRd
installCheckpointRuvPearsonTraces <- function(context, namespace) {
    observer <- function(phase, name, frame) {
        checkpointGuiDiagnosticEvent(context, paste(phase, name, sep = ":"),
            checkpointRuvPearsonShapes(frame))
    }
    context$previous_ruv_observer <- getOption("multischolar.ruv_pearson_observer")
    options(multischolar.ruv_pearson_observer = observer)
    context$ruv_traces <- list()
    add_trace <- function(name, where, signature = NULL, label = name) {
        if (!exists(name, where, inherits = FALSE)) stop(paste("missing RUV/Pearson trace binding:", name))
        enter <- substitute(getOption("multischolar.ruv_pearson_observer")("enter", LABEL, environment()),
            list(LABEL = label))
        leave <- substitute(getOption("multischolar.ruv_pearson_observer")("exit", LABEL, environment()),
            list(LABEL = label))
        arguments <- list(what = name, tracer = enter, exit = leave, print = FALSE, where = where)
        if (!is.null(signature)) arguments$signature <- signature
        # Language-valued tracers must reach trace() unevaluated; otherwise
        # do.call runs the observer during installation, before Y/M exist.
        do.call(base::trace, arguments, quote = TRUE)
        context$ruv_traces[[length(context$ruv_traces) + 1L]] <- list(name = name,
            where = where, signature = signature, label = label)
    }
    for (name in c("resolveProtNormRuvParameters", "findBestNegCtrlPercentage",
        "getRuvIIIReplicateMatrixHelper", "applyProtNormRuvCorrectionStep",
        "finalizeProtNormRuvCleanupStep", "persistProtNormRuvResult",
        "runProtNormCorrelationVectorStep", "runProtNormCorrelationFilterStep")) {
        add_trace(name, namespace)
    }
    for (name in c("runDIAReferenceNormalization", "runDIAReferenceCorrelation")) {
        if (exists(name, namespace, inherits = FALSE)) add_trace(name, namespace)
    }
    for (name in c("getNegCtrlProtAnova", "ruvCancor", "ruvIII_C_Varying", "plotPearson",
        "pearsonCorForSamplePairs", "filterSamplesByProteinCorrelationThreshold")) {
        add_trace(name, namespace, signature = "ProteinQuantitativeData", label = paste0(name, "[ProteinQuantitativeData]"))
    }
    # Tracing an upstream namespace alone can miss its already imported function.
    # These are the actual bindings resolved by the package's frozen methods.
    imports <- parent.env(namespace)
    for (name in c("ruv_cancorplot", "RUVIII_C_Varying", "cor")) {
        add_trace(name, imports, label = paste0("import:", name))
    }
    checkpointGuiDiagnosticEvent(context, "ruv_pearson:traces_installed", list(
        labels = vapply(context$ruv_traces, `[[`, character(1), "label"),
        no_forced_gc = TRUE, no_allocation_profiler = TRUE))
    invisible(NULL)
}

#' Summarize targeted trace coverage without loading the entire diagnostic log.
#' @param directory Completed diagnostic member directory.
#' @return Bounded call counts and first observed matrix shapes; no peak inference.
#' @noRd
readCheckpointRuvPearsonEvidence <- function(directory) {
    connection <- file(file.path(directory, "application-phases.jsonl"), "r")
    on.exit(close(connection), add = TRUE)
    labels <- c("findBestNegCtrlPercentage", "import:RUVIII_C_Varying", "import:ruv_cancorplot",
        "pearsonCorForSamplePairs[ProteinQuantitativeData]", "import:cor",
        "runProtNormCorrelationVectorStep", "runProtNormCorrelationFilterStep")
    owner_path <- file.path(directory, "gui-state-latest.json")
    component_route <- file.exists(owner_path) && file.size(owner_path) <= 65536 &&
        isTRUE(jsonlite::read_json(owner_path)$reference_states)
    if (component_route) labels <- c(setdiff(labels,
        c("runProtNormCorrelationVectorStep", "runProtNormCorrelationFilterStep")),
        "runDIAReferenceNormalization", "runDIAReferenceCorrelation")
    counts <- stats::setNames(lapply(labels, function(x) c(enter = 0L, exit = 0L)), labels)
    shapes <- list()
    complete <- TRUE
    repeat {
        lines <- readLines(connection, n = 128L, warn = FALSE)
        if (!length(lines)) break
        for (line in lines) {
            record <- jsonlite::fromJSON(line, simplifyVector = FALSE)
            if (identical(record$stage, "diagnostic:recording_limit")) complete <- FALSE
            for (phase in c("enter", "exit")) {
                label <- sub(paste0("^", phase, ":"), "", record$stage)
                if (identical(label, record$stage) || !label %in% labels) next
                counts[[label]][[phase]] <- counts[[label]][[phase]] + 1L
                if (phase == "enter" && is.null(shapes[[label]])) shapes[[label]] <- record$detail
            }
        }
    }
    list(coverage_complete = complete, calls = counts, first_call_shapes = shapes,
        timing_eligible = FALSE)
}

#' Start lightweight memory boundaries and mirror stdout without diverting stderr.
#' @param directory Fresh request-owned artifact directory.
#' @return Context without scientific payloads; caller must stop it.
#' @noRd
startCheckpointGuiDiagnostic <- function(directory, collect_heap = FALSE) {
    files <- file.path(directory, c("application-phases.jsonl", "application-console.log"))
    if (!dir.exists(directory) || any(file.exists(files))) stop("existing or invalid diagnostic directory")
    context <- new.env(parent = emptyenv())
    context$directory <- normalizePath(directory)
    context$traced <- character()
    context$sql_queries <- 0
    context$max_log_bytes <- 64 * 1024^2
    context$recording_limited <- FALSE
    context$collect_heap <- isTRUE(collect_heap)
    context$heap_pending <- character()
    context$previous_observer <- getOption("multischolar.gui_diagnostic_observer")
    context$console <- file(files[[2L]], "wt")
    # shinytest2 discovers the server URL through the original stderr stream.
    # Memory boundaries remain durable without tracing every allocation.
    sink(context$console, split = TRUE)
    checkpointGuiDiagnosticEvent(context, "diagnostic:start",
        list(mode = "memory_boundaries", allocation_profile = FALSE, cpu_profile = FALSE,
            max_log_bytes = context$max_log_bytes, completed_action_gc = context$collect_heap,
            timing_eligible = FALSE))
    context
}

#' Stop profiling on normal completion or failure without certifying the journey.
#' @param context Started diagnostic context.
#' @return NULL invisibly.
#' @noRd
stopCheckpointGuiDiagnostic <- function(context) {
    checkpointGuiDiagnosticEvent(context, "diagnostic:stop")
    if (!is.null(context$ruv_traces)) {
        for (entry in rev(context$ruv_traces)) {
            arguments <- list(what = entry$name, where = entry$where)
            if (!is.null(entry$signature)) arguments$signature <- entry$signature
            do.call(base::untrace, arguments)
        }
        options(multischolar.ruv_pearson_observer = context$previous_ruv_observer)
    }
    if (isTRUE(context$sql_traced)) {
        base::untrace("dbSendQuery", signature = c("duckdb_connection", "character"),
            where = asNamespace("duckdb"))
        options(multischolar.gui_sql_observer = context$previous_sql_observer)
    }
    for (name in context$traced) base::untrace(name, where = context$namespace)
    options(multischolar.gui_diagnostic_observer = context$previous_observer)
    sink()
    close(context$console)
    invisible(NULL)
}


#' Verify corresponding repeated outputs without assuming RDS bytes are iteration-invariant.
#' @param directory Completed release-comparison member directory.
#' @param reference Retained single-render reference directory.
#' @return Exact byte-comparison receipt; any missing or mismatched output fails.
#' @noRd
verifyCheckpointAllocatorOutputs <- function(directory, reference) {
    hashes <- lapply(c("control", "trim"), function(arm) lapply(seq_len(2L), function(iteration) {
        paths <- file.path(directory, paste0("allocator-release-cost-", arm),
            paste0("render-", iteration), paste0("proteomics_composite_QC_figure.", c("png", "pdf", "rds")))
        stopifnot(all(file.exists(paths)), all(file.size(paths) > 0))
        stats::setNames(vapply(paths[c(1L, 3L)], function(path)
            digest::digest(file = path, algo = "sha256"), character(1)), c("png", "rds"))
    }))
    expected <- stats::setNames(vapply(file.path(reference,
        paste0("proteomics_composite_QC_figure.", c("png", "rds"))), function(path)
        digest::digest(file = path, algo = "sha256"), character(1)), c("png", "rds"))
    stopifnot(all(vapply(hashes, function(arm) identical(arm[[1L]], expected), logical(1))))
    stopifnot(all(vapply(hashes, function(arm) all(vapply(arm,
        function(value) identical(value[["png"]], expected[["png"]]), logical(1))), logical(1))))
    stopifnot(identical(hashes[[1L]], hashes[[2L]]))
    list(status = "passed_exact_corresponding_outputs", hashes = hashes,
        first_render_matches_reference = TRUE, all_png_match_reference = TRUE,
        matched_iteration_bytes_identical = TRUE,
        iteration_rds_bytes_identical = identical(hashes[[1L]][[1L]][["rds"]],
            hashes[[1L]][[2L]][["rds"]]), performance_acceptance = FALSE)
}
