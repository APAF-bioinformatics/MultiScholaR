# One unchanged mechanism per fresh guarded installed-package process.

#' Fingerprint a completed rollup through bounded independent codec reads.
#' @param owned Existing replay root whose immutable payload is being checked.
#' @param reference Accepted result descriptor.
#' @param samples Optional complete sample IDs for an existing full reference.
#' @return Ordered SHA256 fingerprints, without retaining result batches.
#' @noRd
fingerprintDIAWriterOutput <- function(owned, reference, samples = NULL) {
    identity <- readSavedStoreIdentity(file.path(owned, ".workflow-entities", "store.json"))
    store <- newWorkflowEntityStore(owned, identity$project_id, identity$workflow_id)
    on.exit(closeWorkflowEntityStore(store), add = TRUE)
    # JSON receipts simplify empty lists and descriptor types. The immutable RDS
    # is the authority; verify its identity before opening the ordinary read scope.
    canonical <- readRDS(file.path(store$directory, dirname(reference$payload$relative_path),
        "descriptor.rds"))
    stopifnot(identical(canonical$entity_id, reference$entity_id))
    reference <- canonical
    scope <- newWorkflowTableReadScope(store, reference, reader = "arrow")
    on.exit(releaseWorkflowTableReadScope(scope), add = TRUE, after = FALSE)
    hashes <- character()
    streamSavedWorkflowTable(scope, scope$metadata$columns, function(data, offset, row_names) {
        if (!is.null(samples)) data <- data[as.character(data$Run) %in% samples, , drop = FALSE]
        if (!nrow(data)) return(invisible(NULL))
        attr(data, "dia_source_order") <- NULL
        attr(data, "row.names") <- .set_row_names(nrow(data))
        hashes[[length(hashes) + 1L]] <<- digest::digest(data, algo = "sha256")
    })
    hashes
}

#' Instrument copies of the writer boundary for untimed attribution only.
#' @return Diagnostic bindings and scalar-only observations; no assay owners retained.
#' @noRd
makeDIAWriterDiagnostic <- function() {
    observations <- new.env(parent = emptyenv())
    observations$costs <- list()
    observations$native <- list()
    observations$completed <- observations$failed <- 0L
    observe <- function(label, seconds, bytes = 0, whole = FALSE) {
        previous <- observations$costs[[label]] %||% list(calls = 0L,
            seconds = 0, bytes = 0, whole_batches = 0L)
        observations$costs[[label]] <- list(calls = previous$calls + 1L,
            seconds = previous$seconds + seconds, bytes = previous$bytes + bytes,
            whole_batches = previous$whole_batches + as.integer(whole))
    }
    diagnostic_env <- new.env(parent = environment(appendWorkflowTableWriter))
    diagnostic_env$observe <- observe
    diagnostic_env$observations <- observations
    append <- appendWorkflowTableWriter
    environment(append) <- diagnostic_env
    code <- as.list(body(append))
    matched <- 0L
    for (i in seq_along(code)) {
        expression <- code[[i]]
        if (is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
            is.symbol(expression[[2L]]) &&
            as.character(expression[[2L]]) %in% c("storage", "schema_metadata")) {
            label <- as.character(expression[[2L]])
            code[[i]] <- bquote({
                diagnostic_started <- as.numeric(Sys.time())
                .(expression)
                if (writer$backend == "fast_r") observe(.(label),
                    as.numeric(Sys.time()) - diagnostic_started,
                    as.numeric(object.size(.(as.name(label)))),
                    identical(at, seq_len(nrow(table))))
            })
            matched <- matched + 1L
        }
        if (is.call(expression) && identical(expression[[1L]], as.name("if")) &&
            identical(expression[[2L]], quote(writer$backend == "fast_r"))) {
            branch <- as.list(expression[[3L]])
            branch <- base::append(branch, list(quote(diagnostic_complete <- TRUE)),
                after = length(branch) - 1L)
            expression[[3L]] <- as.call(branch)
            code[[i]] <- expression
        }
    }
    stopifnot(matched == 2L)
    code <- base::append(code, list(quote(diagnostic_complete <- FALSE), quote(on.exit({
        if (writer$backend == "fast_r") {
            if (diagnostic_complete) observations$completed <- observations$completed + 1L else
                observations$failed <- observations$failed + 1L
        }
    }, add = TRUE))), after = 1L)
    body(append) <- as.call(code)
    update <- updateWorkflowFastWriter
    environment(update) <- diagnostic_env
    code <- as.list(body(update))
    code <- base::append(code, list(quote(observations$native[[writer$path]] <- info)), after = 1L)
    body(update) <- as.call(code)
    list(append = append, update = update, observations = observations)
}

#' Make diagnostic-only copies of the query boundary, leaving timed functions intact.
#' @param root Existing guarded member directory.
#' @param observe Scalar boundary observer; never receives assay objects.
#' @return Two bindings for one test scope.
#' @noRd
makeDIARollupDiagnostic <- function(root, observe) {
    original_scope <- newDIARelationalScope
    scope_fn <- function(...) {
        scope <- original_scope(...)
        ready <- FALSE
        on.exit(if (!ready) releaseWorkflowTableReadScope(scope), add = TRUE)
        DBI::dbExecute(scope$connection, "SET profiling_mode = 'detailed'")
        DBI::dbExecute(scope$connection, "SET enable_profiling = 'json'")
        DBI::dbExecute(scope$connection, paste("SET profiling_output =",
            DBI::dbQuoteString(scope$connection, file.path(root, "duckdb-profile.json"))))
        ready <- TRUE
        scope
    }
    stream_fn <- streamDIARelationalQuery
    diagnostic_env <- new.env(parent = environment(stream_fn))
    diagnostic_env$observe <- observe
    diagnostic_env$root <- root
    diagnostic_env$query_number <- 0L
    environment(stream_fn) <- diagnostic_env
    code <- as.list(body(stream_fn))
    matched <- 0L
    for (i in seq_along(code)) {
        expression <- code[[i]]
        if (is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
            identical(expression[[2L]], as.name("result"))) {
            code[[i]] <- quote({
                query_number <<- query_number + 1L
                plan <- DBI::dbGetQuery(con, paste("EXPLAIN", query))
                text <- utils::capture.output(print(plan))
                stopifnot(sum(nchar(text, type = "bytes")) <= 1048576)
                writeLines(text, file.path(root, paste0("query-plan-", query_number, ".txt")))
                observe("before_execute")
                result <- DBI::dbSendQuery(con, query, arrow = TRUE)
                observe("after_execute")
            })
            matched <- matched + 1L
        }
        if (identical(expression, quote(on.exit(DBI::dbClearResult(result), add = TRUE)))) {
            code[[i]] <- quote(on.exit({
                observe("before_clear")
                DBI::dbClearResult(result)
                observe("after_clear")
            }, add = TRUE))
        }
    }
    stopifnot(matched == 1L)
    body(stream_fn) <- as.call(code)
    list(scope = scope_fn, stream = stream_fn)
}

#' Load only the frozen rollup functions, without registering historical S4 methods.
#' @param proof Existing evidence root containing x3Dkg2.
#' @return Frozen helper with its original private function environment.
#' @noRd
bulkRollupFrozenKernel <- function(proof) {
    path <- file.path(proof, "browser-smoke.x3Dkg2", "source", "R", "func_prot_rollup_methods.R")
    frozen <- new.env(parent = asNamespace("MultiScholaR"))
    helpers <- c(".rollupSumObservedOrNA", ".rollupAssertPopulatedKey", ".rollupAssertLinearQuantity",
        ".rollupDuplicateExamples", ".rollupSummarisePrepared", "rollUpPrecursorToPeptideHelper")
    for (expr in parse(path)) {
        if (is.call(expr) && identical(expr[[1L]], as.name("<-")) && is.symbol(expr[[2L]]) &&
            as.character(expr[[2L]]) %in% helpers) eval(expr, frozen)
    }
    stopifnot(all(vapply(helpers, exists, logical(1), envir = frozen, inherits = FALSE)))
    frozen$rollUpPrecursorToPeptideHelper
}

#' Assemble the experimental bulk input through the existing bounded reader/buffer.
#' @param scope Validated source scope.
#' @param input Declared rollup roles and optional evidence/annotation.
#' @param samples Complete sample identities, or NULL for the full population.
#' @return One table, in original source order, with no retained batch list.
#' @noRd
readDIABulkRollup <- function(scope, input, samples = NULL) {
    fields <- unique(c(unlist(input$roles), input$evidence, input$annotation))
    selected <- NULL
    if (!is.null(samples)) {
        index <- saveWorkflowGroupIndex(scope, input$roles$sample)
        groups <- readWorkflowIndexDirectory(scope$store, index$value$groups)
        groups <- groups[match(samples, as.character(groups$value)), , drop = FALSE]
        stopifnot(!anyNA(groups$size))
        index_scope <- newWorkflowTableReadScope(scope$store, index$value$positions, reader = "arrow")
        on.exit(releaseWorkflowTableReadScope(index_scope), add = TRUE)
        index_rows <- sort(unlist(Map(seq.int, groups$start, length.out = groups$size)))
        selected <- integer(sum(groups$size))
        streamSavedWorkflowTable(index_scope, "source_order", function(data, offset, row_names) {
            selected[offset + seq_len(nrow(data))] <<- data$source_order
        }, source_rows = index_rows)
        releaseWorkflowTableReadScope(index_scope)
        selected <- sort(selected)
    }
    rows <- if (is.null(samples)) scope$metadata$rows else length(selected)
    buffer <- newWorkflowTableBuffer(scope$metadata, fields, rows)
    streamSavedWorkflowTable(scope, fields, function(data, offset, row_names) {
        at <- offset + seq_len(nrow(data))
        source_order <- if (is.null(selected)) at else attr(data, "dia_source_order")
        assignWorkflowTableBuffer(buffer, data, at, source_order, row_names)
    }, source_rows = selected)
    table <- finishWorkflowTableBuffer(buffer)
    attr(table, "workflow_source_order") <- NULL
    table
}

#' Execute and persist one native bulk call using unchanged scientific kernels.
#' @param store,input,persist Existing rollup adapter arguments.
#' @param kernel Frozen or current helper; never traced or wrapped internally.
#' @param cluster Eight baseline workers or NA for the current serial branch.
#' @param samples Complete sample selector or NULL.
#' @param phase Natural-boundary timing callback.
#' @return Existing reference, summary, dependency and filtering-count contract.
#' @noRd
computeDIABulkRollup <- function(store, input, persist, kernel, cluster, samples, phase) {
    scope <- newWorkflowTableReadScope(store, input$reference, reader = "arrow")
    on.exit(releaseWorkflowTableReadScope(scope), add = TRUE)
    table <- phase("read", function() readDIABulkRollup(scope, input, samples))
    releaseWorkflowTableReadScope(scope)
    roles <- input$roles
    result <- phase("kernel", function() kernel(table,
        sample_id_column = roles$sample, protein_id_column = roles$protein,
        peptide_sequence_column = roles$peptide, modified_peptide_sequence_column = roles$modified,
        precursor_quantity_column = roles$raw, precursor_normalised_column = roles$normalised,
        precursor_id_column = roles$precursor, precursor_charge_column = "Precursor.Charge",
        core_utilisation = cluster, is_logged_data = FALSE, return_rollup_result = TRUE))
    table <- NULL
    counts <- newDIAFilteringAccumulator(persist$prototype, roles$protein)
    metadata <- workflowSavedTableMetadata(result$data[0, , drop = FALSE], ".rownames", 3L)
    metadata$rows <- nrow(result$data)
    reference <- phase("summary_and_write", function() saveWorkflowTable(store, NULL,
        basename(tempfile("bulk-rollup-")), persist$identity, metadata$key_column,
        metadata = metadata, write_fn = function(path, restoration) {
            writer <- openWorkflowTableWriter(path, restoration, store$progress)
            on.exit(abortWorkflowTableWriter(writer), add = TRUE)
            for (start in seq.int(1L, max(1L, nrow(result$data)), by = 10000L)) {
                at <- if (nrow(result$data)) seq.int(start, min(nrow(result$data), start + 9999L)) else integer()
                counts$push(result$data[at, , drop = FALSE])
                appendWorkflowTableWriter(writer, result$data, at)
            }
            finishWorkflowTableWriter(writer)
        }))
    list(reference = reference, summary = result$summary, columns = names(result$data),
        dependencies = list(input$reference), counts = counts$finish())
}

#' Run one complete bulk operation with the existing checkpoint and publication owners.
#' @param manager Owned replay adapter positioned at protein_peptide_filtered.
#' @param root,proof Member output and retained evidence roots.
#' @param arm A (frozen eight workers) or B (current single process).
#' @param samples Complete sample selector or NULL.
#' @param phase Phase timing callback.
#' @param failure Inject publication failure for the small lifecycle proof.
#' @return Accepted descriptor and required-output inventory.
#' @noRd
runDIABulkOperation <- function(manager, root, proof, arm, samples, phase, failure = FALSE) {
    profile <- if (grepl("^P([0-9]|1[0-6])$", arm)) arm else "P0"
    if (grepl("^P([0-9]|1[0-6])$", arm)) arm <- "D"
    if (identical(Sys.getenv("MULTISCHOLAR_PERF108_LOCAL"), "true") && profile == "P14") arm <- "C"
    if (identical(Sys.getenv("MULTISCHOLAR_PERF108_NATIVE"), "true") && profile == "P14") arm <- "F"
    cluster <- NA
    on.exit(if (inherits(cluster, "multidplyr_cluster")) stopDIAWorkerCluster(cluster), add = TRUE)
    kernel <- if (arm == "A") bulkRollupFrozenKernel(proof) else rollUpPrecursorToPeptideHelper
    if (arm == "A") cluster <- phase("worker_start", function() startDIAWorkerCluster(newDIAWorkerSpec(8L)))
    workers <- if (arm == "A") vapply(cluster, \(worker) worker$get_pid(), integer(1)) else Sys.getpid()
    testthat::local_mocked_bindings(computeSavedWorkflowPrecursorRollup = function(store, input, persist = NULL) {
        stopifnot(!is.null(persist))
        if (arm %in% c("C", "D", "F")) return(computeDIARelationalRollup(store, input, persist,
            arm, samples, phase, if (arm == "F") NULL else diaRollupExtensionInfo()$path, profile = profile))
        computeDIABulkRollup(store, input, persist, kernel, cluster, samples, phase)
    })
    project <- file.path(root, "dia_baseline_smoke")
    dir.create(project)
    directories <- setupDirectories(project, "proteomics")
    directories$proteomics_dia_baseline_smoke <- directories$proteomics
    manager$.context$experiment_paths <- directories$proteomics
    manager$.context$reference_states <- TRUE
    manager$.context$session_format <- "native_v1"
    config <- readDIAStateMetadata(manager)$config
    rlang::local_bindings(config_list = config, project_dirs = directories,
        filtering_progress = methods::new("FilteringProgress"), .env = .GlobalEnv)
    withr::local_options(multischolar.test_mode = FALSE)
    withr::local_envvar(MULTISCHOLAR_TEST_MODE = "false")
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- manager
    workflow$qc_params <- list(peptide_qc = list())
    grDevices::pdf(NULL)
    device <- grDevices::dev.cur()
    on.exit(if (device %in% grDevices::dev.list()) grDevices::dev.off(device), add = TRUE)
    result <- withDIAStateScope(manager, function() {
        result <- phase("apply_including_read_kernel_write_commit", function() runDIAReferencePrecursorRollup(workflow))
        arguments <- list(data = result$filteringSummary, step_name = "4_precursor_rollup",
            omic_type = "proteomics", experiment_label = "dia_baseline_smoke",
            return_grid = TRUE, overwrite = TRUE, project_dirs = directories)
        if (failure) {
            expect_error(phase("injected_publication_failure", function()
                publishDIAQvalueOutput(manager, result, arguments,
                    function(...) stop("injected bulk publication failure"))), "injected bulk publication failure")
        }
        phase("publication", function() publishDIAQvalueOutput(manager, result, arguments,
            updateProteinFiltering))
        result
    })
    if (arm == "A") {
        phase("worker_stop", function() stopDIAWorkerCluster(cluster))
        stopifnot(all(vapply(cluster, \(worker) !worker$is_alive(), logical(1))))
        cluster <- NA
    }
    record <- readDIAStateMetadata(manager)
    files <- list.files(directories$proteomics$time_dir, recursive = TRUE, full.names = TRUE)
    figures <- files[grepl("4_precursor_rollup_.*[.](png|pdf)$", files) & !grepl("/figures/", files)]
    stopifnot(length(figures) == 12L, all(file.size(figures) > 0), manager$.context$store$open_resources == 0L)
    list(reference = record$data$components$peptide_data,
        summary = record$data$slots$args$rollUpPrecursorToPeptide$rollup_summary,
        filtering = unclass(result$filteringSummary), worker_pids = workers,
        figures = stats::setNames(vapply(figures, file.size, numeric(1)), basename(figures)),
        state = result$stateReference$entity_id, execution = manager$.context$last_rollup_execution)
}

#' Create an isolated replay using immutable payload links and independent metadata.
#' @param source_store Frozen committed store directory.
#' @param destination Fresh request-owned directory on the same filesystem.
#' @return Reopened reference adapter; original pending files are excluded.
#' @noRd
copyDIAReplayCheckpoint <- function(source_store, destination) {
    stopifnot(dir.exists(source_store), dir.exists(destination),
        !length(list.files(destination, all.files = TRUE, no.. = TRUE)))
    target <- file.path(destination, ".workflow-entities")
    stopifnot(dir.create(target))
    for (name in c("store.json", "dia-state-head.rds"))
        stopifnot(file.copy(file.path(source_store, name), target))
    files <- list.files(file.path(source_store, "entities"), recursive = TRUE, full.names = TRUE)
    files <- files[!dir.exists(files)]
    for (path in files) {
        relative <- substring(path, nchar(source_store) + 2L)
        to <- file.path(target, relative)
        dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
        if (grepl("/payload[.](parquet|rds|qs2)$", path)) stopifnot(file.link(path, to)) else
            stopifnot(file.copy(path, to))
    }
    identity <- readSavedStoreIdentity(file.path(target, "store.json"))
    newDIAStateAdapter(newWorkflowEntityStore(destination, identity$project_id, identity$workflow_id), reopen = TRUE)
}


#' Rebuild and publish progress for a selected qualification population.
#' @param manager,root Owned selected checkpoint and output directory.
#' @param failure Exercise failed publication and retry without changing science.
#' @return Published progress bound to the selected table and design.
#' @noRd
publishDIARollupQualification <- function(manager, root, failure = FALSE) {
    project <- file.path(root, "qualification-publication")
    dir.create(project)
    directories <- setupDirectories(project, "proteomics")
    directories$proteomics_dia_baseline_smoke <- directories$proteomics
    manager$.context$experiment_paths <- directories$proteomics
    manager$.context$reference_states <- TRUE
    record <- readDIAStateMetadata(manager)
    rlang::local_bindings(config_list = record$config, project_dirs = directories,
        filtering_progress = methods::new("FilteringProgress"), .env = .GlobalEnv)
    withr::local_options(multischolar.test_mode = FALSE)
    withr::local_envvar(MULTISCHOLAR_TEST_MODE = "false")
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- manager
    workflow$qc_params <- list(peptide_qc = list())
    grDevices::pdf(NULL)
    device <- grDevices::dev.cur()
    on.exit(if (device %in% grDevices::dev.list()) grDevices::dev.off(device), add = TRUE)
    withDIAStateScope(manager, function() {
        expect_error(assertDIAExportReady(manager), class = "saved_entity_outputs_pending")
        parameters <- record$data$slots$args$filterMinNumPeptidesPerProtein[
            c("num_peptides_per_protein_thresh", "num_peptidoforms_per_protein_thresh")]
        stopifnot(all(vapply(parameters, function(value) is.numeric(value) &&
            length(value) == 1L && is.finite(value), logical(1))))
        result <- finishDIAReferencePeptideStep(workflow, record$data,
            "protein_peptide_filtered", "selected qualification population", parameters,
            "filterMinNumPeptidesPerProtein")
        expect_error(assertDIAExportReady(manager), class = "saved_entity_outputs_pending")
        arguments <- list(data = result$filteringSummary, step_name = "3_protein_peptide_filter",
            omic_type = "proteomics", experiment_label = "dia_baseline_smoke",
            return_grid = TRUE, overwrite = TRUE, project_dirs = directories)
        if (failure) {
            expect_error(publishDIAQvalueOutput(manager, result, arguments,
                function(...) stop("injected qualification publication failure")),
                "injected qualification publication failure")
            expect_error(assertDIAExportReady(manager), class = "saved_entity_outputs_pending")
        }
        publishDIAQvalueOutput(manager, result, arguments, updateProteinFiltering)
        expect_null(assertDIAExportReady(manager))
        result
    })
}

#' Repair only affected prepared checkpoints without regenerating valid inputs.
#' @param proof,root Existing evidence root and fresh guarded output directory.
#' @param preflight Repair the smallest affected case and prove actual Apply.
#' @return Revision-bound manifest; original artifacts remain immutable.
#' @noRd
repairDIARollupQualification <- function(proof, root, preflight = FALSE) {
    original <- file.path(proof, "browser-smoke.Ny7QzY", "rollup-qualification-inputs.json")
    manifest <- jsonlite::read_json(original)
    stopifnot(identical(manifest$status, "complete"))
    manifest$revision <- "p3-adaptations-v1"
    manifest$predecessor_sha256 <- digest::digest(file = original, algo = "sha256")
    manifest$status <- "running"
    path <- file.path(root, "rollup-qualification-inputs.json")
    flush <- function() jsonlite::write_json(manifest, path, auto_unbox = TRUE,
        pretty = TRUE, null = "null", digits = NA)
    on.exit(flush(), add = TRUE)
    tryCatch({
        for (case in if (preflight) "9600" else c("9600", "96000", "960000")) {
            selected <- manifest$cases[[case]]
            head <- file.path(selected$store_directory, "dia-state-head.rds")
            stopifnot(identical(digest::digest(file = head, algo = "sha256"), selected$head_sha256))
            owned <- file.path(root, paste0("size-", case))
            dir.create(owned)
            local({
                manager <- copyDIAReplayCheckpoint(selected$store_directory, owned)
                on.exit(closeDIAStateAdapter(manager), add = TRUE)
                manager$.context$reference_states <- TRUE
                result <- publishDIARollupQualification(manager, owned, failure = preflight)
                expect_identical(readDIAStateMetadata(manager)$data$components$peptide_data$entity_id,
                    selected$reference$entity_id)
                expect_identical(manager$.context$store$open_resources, 0L)
                if (preflight) {
                    value <- runDIABulkOperation(manager, owned, proof, "P3", NULL,
                        function(label, fn) fn(), failure = TRUE)
                    expect_identical(value$summary$input_precursor_rows, 9600L)
                    expect_identical(value$summary$output_stripped_peptide_identities, 4800L)
                }
                updated <- selected
                updated$root <- owned
                updated$store_directory <- manager$.context$store$directory
                updated$head_sha256 <- digest::digest(file = file.path(updated$store_directory,
                    "dia-state-head.rds"), algo = "sha256")
                manifest$cases[[case]] <<- updated
            })
            stopifnot(identical(digest::digest(file = head, algo = "sha256"), selected$head_sha256))
            flush()
        }
        manifest$status <- if (preflight) "preflight_complete" else "complete"
    }, error = function(error) {
        manifest$status <<- "failed"
        manifest$error <<- conditionMessage(error)
        stop(error)
    })
    flush()
    manifest
}

#' Prepare immutable rollup qualification inputs with the existing store and writer.
#' @param proof,root Existing evidence root and fresh guarded member directory.
#' @return Manifest of case roots, input hashes and bounded preparation timings.
#' @noRd
prepareDIARollupQualification <- function(proof, root, preflight = FALSE, inject_failure = FALSE) {
    stopifnot(!inject_failure || preflight)
    total_proteins <- if (preflight) 8L else 10000L
    contrast_rows <- total_proteins * 384L
    scope <- manager <- NULL
    completed <- FALSE
    source <- jsonlite::read_json(file.path(proof, "browser-smoke.BpSyxB", "gui-state-latest.json"))
    source_hash <- digest::digest(file = file.path(source$store_directory, "dia-state-head.rds"), algo = "sha256")
    stopifnot(identical(source_hash, "5d64ca562669d94eea1677a1a03d5b410a8a81d8ea59add8cabd6065deb91353"))
    identity <- list(input = "rollup-qualification", design = "frozen", parameters = "20260929",
        selectors = "declared-complete-groups", software = "DNMnfd")
    manifest <- list(status = "running", source_head_sha256 = source_hash, cases = list())
    on.exit({
        if (!is.null(scope) && !isTRUE(scope$released)) releaseWorkflowTableReadScope(scope)
        if (!is.null(manager)) closeDIAStateAdapter(manager)
        if (!completed) manifest$status <- "failed"
        jsonlite::write_json(manifest, file.path(root, "rollup-qualification-inputs.json"),
            auto_unbox = TRUE, pretty = TRUE, null = "null")
    }, add = TRUE)
    record <- function(manager, case, owned, started) {
        spec <- readDIAStateMetadata(manager)$data
        head <- file.path(manager$.context$store$directory, "dia-state-head.rds")
        value <- list(root = owned, store_directory = manager$.context$store$directory,
            rows = spec$dimensions$rows, reference = spec$components$peptide_data,
            head_sha256 = digest::digest(file = head, algo = "sha256"),
            preparation_seconds = proc.time()[["elapsed"]] - started)
        jsonlite::write_json(list(store_directory = value$store_directory),
            file.path(owned, "gui-state-latest.json"), auto_unbox = TRUE)
        manifest$cases[[case]] <<- value
        closeDIAStateAdapter(manager)
        jsonlite::write_json(manifest, file.path(root, "rollup-qualification-inputs.json"),
            auto_unbox = TRUE, pretty = TRUE, null = "null")
    }
    for (rows in if (preflight) integer() else c(4608L, 9600L, 96000L, 960000L)) {
        started <- proc.time()[["elapsed"]]
        owned <- file.path(root, paste0("size-", rows))
        dir.create(owned)
        if (rows == 4608L) {
            manager <- newDIAStateAdapter(newWorkflowEntityStore(owned, "crossover-small", "DIA"),
                WorkflowState$new(audit_enabled = FALSE))
            manager$.context$reference_states <- TRUE
            manager$saveState("protein_peptide_filtered",
                makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small")), list(), "small crossover")
        } else {
            manager <- copyDIAReplayCheckpoint(source$store_directory, owned)
            selectDIAStateReference(manager, "protein_peptide_filtered")
            before <- readDIAStateMetadata(manager)
            checkpoint <- before$data
            scope <- newWorkflowTableReadScope(manager$.context$store, checkpoint$components$peptide_data)
            con <- scope$connection
            DBI::dbExecute(con, "SET threads = 8")
            DBI::dbExecute(con, "SET memory_limit = '12GiB'")
            design <- withDIAStateScope(manager, function()
                readDIAExecutionComponents(manager, fields = "design_matrix"))@design_matrix
            samples <- as.character(design$Run[c(1:4, 149:152, 297:300)])
            proteins <- as.integer(rows / 96L)
            q <- function(x) as.character(DBI::dbQuoteIdentifier(con, x))
            fields <- workflowTableStorageColumns(scope$metadata, scope$metadata$columns)
            input_path <- DBI::dbQuoteString(con, scope$path)
            predicate <- paste0(q("Run"), " IN (", paste(DBI::dbQuoteString(con, samples), collapse = ","),
                ") AND ", q("Protein.Group"), " IN (SELECT DISTINCT ", q("Protein.Group"),
                " FROM read_parquet(", input_path, ") ORDER BY ", q("Protein.Group"), " LIMIT ", proteins, ")")
            generated_key <- if (scope$metadata$key_column == ".rownames" && isTRUE(scope$metadata$automatic_row_names))
                ', CAST(row_number() OVER (ORDER BY ".multischolar_order") AS VARCHAR) AS ".multischolar_row_key"' else ""
            query <- paste0("SELECT ", paste(q(fields), collapse = ","), generated_key,
                ', CAST(row_number() OVER (ORDER BY ".multischolar_order") AS INTEGER) AS ".multischolar_order"',
                " FROM read_parquet(", input_path, ") WHERE ", predicate, ' ORDER BY ".multischolar_order"')
            metadata <- scope$metadata
            metadata$rows <- rows
            selected <- saveWorkflowTableQuery(scope, query, metadata, paste0("selection-", rows), identity)
            releaseWorkflowTableReadScope(scope)
            checkpoint$components$peptide_data <- selected
            checkpoint$dimensions <- list(rows = rows, columns = metadata$columns)
            withDIAStateScope(manager, function() {
                checkpoint$components$design_matrix <- saveDIATableComponent(manager$.context,
                    design[as.character(design$Run) %in% samples, , drop = FALSE])
                commitDIAReferenceState(manager, "protein_peptide_filtered", checkpoint,
                    list(qualification_rows = rows), "complete-group crossover selection")
            })
            publishDIARollupQualification(manager, owned)
        }
        record(manager, as.character(rows), owned, started)
    }
    started <- proc.time()[["elapsed"]]
    owned <- file.path(root, "contrast")
    dir.create(owned)
    manager <- newDIAStateAdapter(newWorkflowEntityStore(owned, "rollup-contrast", "DIA"),
        WorkflowState$new(audit_enabled = FALSE))
    manager$.context$reference_states <- TRUE
    recipe <- readWorkflowCheckpointRecipe("annotation_heavy")
    recipe$peptides_per_protein <- 2L
    make_chunk <- function(proteins) {
        table <- makeWorkflowPrecursorChunk(recipe, proteins)
        table <- table[rep(seq_len(nrow(table)), each = 2L), , drop = FALSE]
        table$Precursor.Charge <- rep(c(2L, 3L), nrow(table) / 2L)
        table$Precursor.Id <- paste(table$Modified.Sequence, table$Precursor.Charge, sep = "_")
        labels <- vapply(proteins, function(protein) paste0(
            paste(sprintf("A%06X", (protein - 1L) * 256L + 0:255), collapse = ";"), ";"), character(1))
        stopifnot(all(nchar(labels, type = "bytes") == 2048L))
        table$Protein.Ids <- labels[match(table$Protein.Group, sprintf("PG_%06d", proteins))]
        table$identification_peptide_count <- table$peptides_for_protein_count <- 2L
        table$identification_peptidoform_count <- table$peptidoforms_for_protein_count <- 4L
        rownames(table) <- NULL
        table
    }
    prototype <- make_chunk(1L)[0, , drop = FALSE]
    metadata <- workflowSavedTableMetadata(prototype, ".rownames", 3L)
    metadata$rows <- contrast_rows
    reference <- saveWorkflowTable(manager$.context$store, NULL, "contrast-input", identity,
        metadata$key_column, metadata = metadata, write_fn = function(path, restoration) {
            writer <- openWorkflowTableWriter(path, restoration)
            on.exit(abortWorkflowTableWriter(writer), add = TRUE)
            for (start in seq.int(1L, total_proteins, by = 256L)) {
                chunk <- make_chunk(seq.int(start, min(total_proteins, start + 255L)))
                if (inject_failure) stop("injected qualification preparation failure")
                for (offset in seq.int(1L, nrow(chunk), by = 10000L))
                    appendWorkflowTableWriter(writer, chunk, seq.int(offset, min(nrow(chunk), offset + 9999L)))
            }
            finishWorkflowTableWriter(writer)
        })
    object <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    object@peptide_data <- prototype
    object@design_matrix <- makeWorkflowCheckpointDesign(48L)
    object@args$fixture <- list(recipe = "rollup-contrast-20260929", seed = 902005L, synthetic = TRUE,
        n_proteins = total_proteins, n_samples = 48L, peptides_per_protein = 2L, modified_forms = 2L,
        charges = 2L, annotation_bytes = 2048L)
    manager$saveState("protein_peptide_filtered", object, list(), "contrast prototype")
    checkpoint <- readDIAStateMetadata(manager)$data
    checkpoint$components$peptide_data <- reference
    checkpoint$dimensions <- list(rows = contrast_rows, columns = names(prototype))
    withDIAStateScope(manager, function() commitDIAReferenceState(manager, "protein_peptide_filtered",
        checkpoint, list(seed = 902005L, forms = 2L, charges = 2L, annotation_bytes = 2048L),
        "frozen contrasting rollup input"))
    record(manager, "contrast", owned, started)
    manifest$status <- "complete"
    completed <- TRUE
    stopifnot(identical(digest::digest(file = file.path(source$store_directory, "dia-state-head.rds"),
        algo = "sha256"), source_hash))
    jsonlite::write_json(manifest, file.path(root, "rollup-qualification-inputs.json"),
        auto_unbox = TRUE, pretty = TRUE, null = "null")
    manifest
}

#' Prepare existing wide-projection and ordinal-index I/O workloads.
#' @param proof,root Existing evidence root and fresh guarded directory.
#' @return Frozen references and native object fingerprints for bounded requests.
#' @noRd
prepareDIARollupIO <- function(proof, root) {
    owned <- file.path(root, "store")
    dir.create(owned)
    source <- jsonlite::read_json(file.path(proof, "browser-smoke.BpSyxB", "gui-state-latest.json"))
    manager <- copyDIAReplayCheckpoint(source$store_directory, owned)
    on.exit(closeDIAStateAdapter(manager), add = TRUE)
    input <- readDIAStateMetadata(manager, getDIAStateReference(manager, "protein_peptide_filtered"))$data$components$peptide_data
    scope <- newWorkflowTableReadScope(manager$.context$store, input, reader = "arrow")
    on.exit(releaseWorkflowTableReadScope(scope), add = TRUE, after = FALSE)
    index <- saveWorkflowGroupIndex(scope, "Run")$value$positions
    releaseWorkflowTableReadScope(scope)
    protein <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("wide"), allow_large = TRUE,
        output_root = root, repo_root = here::here())
    table <- protein@protein_quant_table
    identity <- list(input = "wide-902004", design = "256-samples", parameters = "bounded-io",
        selectors = "all", software = "bvT5nP")
    wide <- saveWorkflowTable(manager$.context$store, table, "wide-io", identity,
        ".rownames", metadata_version = 3L)
    requests <- list(index = list(offset = c(0L, 11995000L, 23990000L), limit = c(10000L, 10000L, 10000L)),
        wide = list(offset = c(0L, 3500L, 7000L), limit = c(25L, 1000L, 1000L)))
    expected_wide <- vapply(seq_along(requests$wide$offset), function(i) {
        rows <- requests$wide$offset[[i]] + seq_len(requests$wide$limit[[i]])
        fingerprintDIARollupIO(restoreWorkflowAttributeOrder(table[rows, , drop = FALSE],
            names(attributes(table))))
    }, character(1))
    # Independent Arrow/buffer oracle for exactly the bounded index requests.
    index_scope <- newWorkflowTableReadScope(manager$.context$store, index, reader = "arrow")
    on.exit(releaseWorkflowTableReadScope(index_scope), add = TRUE, after = FALSE)
    expected_index <- vapply(seq_along(requests$index$offset), function(i) {
        rows <- requests$index$offset[[i]] + seq_len(requests$index$limit[[i]])
        buffer <- newWorkflowTableBuffer(index_scope$metadata, "source_order", length(rows))
        streamSavedWorkflowTable(index_scope, "source_order", function(data, offset, row_names) {
            assignWorkflowTableBuffer(buffer, data, offset + seq_len(nrow(data)),
                attr(data, "dia_source_order"), row_names)
        }, source_rows = rows)
        data <- finishWorkflowTableBuffer(buffer)
        attr(data, "workflow_source_order") <- NULL
        fingerprintDIARollupIO(data)
    }, character(1))
    releaseWorkflowTableReadScope(index_scope)
    value <- list(status = "complete", root = owned, identity = manager$.context$store$identity,
        index = index, wide = wide, requests = requests,
        expected = list(index = expected_index, wide = expected_wide))
    closeDIAStateAdapter(manager)
    saveRDS(value, file.path(root, "rollup-io-inputs.rds"), compress = FALSE)
    value
}

#' Correct the I/O oracle's table attribute order while retaining its input files.
#' @param proof,root Existing evidence root and a fresh guarded preflight directory.
#' @return New oracle revision referencing the unchanged prepared Parquet artifacts.
#' @noRd
repairDIARollupIOInputs <- function(proof, root) {
    path <- file.path(proof, "browser-smoke.QnPzsh", "rollup-io-inputs.rds")
    input <- readRDS(path)
    stopifnot(identical(input$status, "complete"))
    protein <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("wide"), allow_large = TRUE,
        output_root = root, repo_root = here::here())
    table <- protein@protein_quant_table
    store <- newWorkflowEntityStore(input$root, input$identity$project_id, input$identity$workflow_id)
    on.exit(closeWorkflowEntityStore(store), add = TRUE)
    scope <- newWorkflowTableReadScope(store, input$wide, reader = "arrow")
    on.exit(releaseWorkflowTableReadScope(scope), add = TRUE, after = FALSE)
    stopifnot(identical(scope$metadata$columns, names(table)),
        identical(scope$metadata$attribute_order, names(attributes(table))))
    input$expected$wide <- vapply(seq_along(input$requests$wide$offset), function(i) {
        rows <- input$requests$wide$offset[[i]] + seq_len(input$requests$wide$limit[[i]])
        fingerprintDIARollupIO(restoreWorkflowAttributeOrder(table[rows, , drop = FALSE],
            scope$metadata$attribute_order))
    }, character(1))
    input$revision <- "codec-attributes-v2"
    input$predecessor_sha256 <- digest::digest(file = path, algo = "sha256")
    input
}

#' Exercise the existing bounded access contract with an explicit eight-thread reader.
#' @param input Frozen I/O preparation record.
#' @param workload,reader Either index/wide and duckdb/fastduckdb.
#' @param root Existing guarded output directory.
#' @return A small measurement and exactness receipt.
#' @noRd
runDIARollupIO <- function(input, workload, reader, root, failure = FALSE) {
    started <- proc.time()
    store <- newWorkflowEntityStore(input$root, input$identity$project_id, input$identity$workflow_id)
    on.exit(closeWorkflowEntityStore(store), add = TRUE)
    scope <- newWorkflowTableReadScope(store, input[[workload]], reader = reader)
    on.exit(releaseWorkflowTableReadScope(scope), add = TRUE, after = FALSE)
    if (reader == "duckdb") {
        DBI::dbExecute(scope$connection, "SET threads = 8")
        DBI::dbExecute(scope$connection, "SET memory_limit = '12GiB'")
    }
    read_fn <- readSavedWorkflowTable
    if (reader == "fastduckdb") {
        code <- deparse(body(read_fn), width.cutoff = 500L)
        stopifnot(sum(grepl("n_threads = 2L", code, fixed = TRUE)) == 1L)
        body(read_fn) <- parse(text = gsub("n_threads = 2L", "n_threads = 8L", code, fixed = TRUE))[[1L]]
    }
    columns <- scope$metadata$columns
    blocks <- split(columns, ceiling(seq_along(columns) / 64L))
    on.exit({
        releaseWorkflowTableReadScope(scope)
        stopifnot(store$open_resources == 0L)
        jsonlite::write_json(list(reader = reader, failure = failure, open_resources = store$open_resources),
            file.path(root, paste0("io-cleanup-", reader, "-", workload, "-", failure, ".json")), auto_unbox = TRUE)
    }, add = TRUE, after = FALSE)
    results <- character(length(input$requests[[workload]]$offset))
    for (i in seq_along(results)) {
        values <- list()
        for (block in blocks) {
            page <- read_fn(store, input[[workload]], block,
                offset = input$requests[[workload]]$offset[[i]],
                limit = input$requests[[workload]]$limit[[i]], scope = scope, reader = reader)
            stopifnot(identical(rownames(page$data), as.character(
                input$requests[[workload]]$offset[[i]] + seq_len(nrow(page$data)))))
            if (failure) stop("injected bounded reader failure")
            values[block] <- unclass(page$data)
        }
        for (attribute in names(scope$metadata$attributes))
            attr(values, attribute) <- scope$metadata$attributes[[attribute]]
        attr(values, "row.names") <- .set_row_names(nrow(page$data))
        class(values) <- scope$metadata$class
        values <- restoreWorkflowAttributeOrder(values, scope$metadata$attribute_order)
        results[[i]] <- fingerprintDIARollupIO(values)
    }
    stopifnot(identical(unname(results), unname(input$expected[[workload]])))
    releaseWorkflowTableReadScope(scope)
    closeWorkflowEntityStore(store)
    list(status = "complete", reader = reader, threads = 8L, workload = workload,
        seconds = unname((proc.time() - started)[["elapsed"]]),
        cpu = unname(sum((proc.time() - started)[c("user.self", "sys.self")])),
        hashes = results, reference = input[[workload]]$entity_id,
        scope = "bounded projection/range consumer; not rollup computation")
}

#' Fingerprint bounded I/O values, preserving finite bits and NA/NaN distinctions.
#' @param data Restored bounded table; never an assay-sized input.
#' @return SHA256 with only automatic block row names and non-NaN NA payloads canonicalized.
#' @noRd
fingerprintDIARollupIO <- function(data) {
    attr(data, "row.names") <- .set_row_names(nrow(data))
    for (name in names(data)) if (is.double(data[[name]])) {
        # R arithmetic can quiet the NA payload; Parquet restores canonical NA.
        # This is an identical R missing value, not a numerical tolerance.
        missing <- is.na(data[[name]]) & !is.nan(data[[name]])
        if (any(missing)) data[[name]][missing] <- NA_real_
    }
    digest::digest(data, algo = "sha256")
}

#' Inspect physical sample pruning without executing a rollup or altering inputs.
#' @param proof,root Existing evidence and a fresh guarded output directory.
#' @return Bounded source statistics and aggregation-only physical plans.
#' @noRd
inspectDIARollupScans <- function(proof, root) {
    source <- jsonlite::read_json(file.path(proof, "browser-smoke.BpSyxB", "gui-state-latest.json"))
    original_head <- file.path(source$store_directory, "dia-state-head.rds")
    original_hash <- digest::digest(file = original_head, algo = "sha256")
    owned <- file.path(root, "inspection")
    dir.create(owned)
    manager <- copyDIAReplayCheckpoint(source$store_directory, owned)
    on.exit(closeDIAStateAdapter(manager), add = TRUE)
    selectDIAStateReference(manager, "protein_peptide_filtered")
    checkpoint <- readDIAStateMetadata(manager)$data
    input <- list(reference = checkpoint$components$peptide_data,
        roles = list(sample = "Run", protein = "Protein.Group", peptide = "Stripped.Sequence",
            modified = "Modified.Sequence", precursor = "Precursor.Id",
            raw = "Precursor.Quantity", normalised = "Precursor.Normalised"))
    scope <- newDIARelationalScope(manager$.context$store, input, profile = "P3")
    on.exit(releaseWorkflowTableReadScope(scope), add = TRUE, after = FALSE)
    con <- scope$connection
    q <- function(x) as.character(DBI::dbQuoteIdentifier(con, x))
    DBI::dbExecute(con, paste0("CREATE TEMP VIEW qc_input AS SELECT ",
        paste(q(c(unlist(input$roles), ".multischolar_order")), collapse = ","),
        " FROM read_parquet(", DBI::dbQuoteString(con, scope$path), ")"))
    samples <- DBI::dbGetQuery(con, 'SELECT DISTINCT "Run" FROM qc_input ORDER BY "Run"')$Run
    stopifnot(length(samples) == 300L)
    plans <- list()
    for (profile in c("P0", "P6", "P3")) {
        selected <- if (profile == "P3") samples[1:100] else samples
        predicate <- if (profile == "P6") "" else paste0(' WHERE "Run" IN (',
            paste(DBI::dbQuoteString(con, selected), collapse = ","), ")")
        DBI::dbExecute(con, paste0("CREATE OR REPLACE TEMP VIEW qc_aggregate_input AS ",
            "SELECT * FROM qc_input", predicate))
        prepareDIARelationalAggregate(scope, input, "D")
        plans[[profile]] <- DBI::dbGetQuery(con,
            'EXPLAIN SELECT * FROM qc_aggregated ORDER BY "Run", "Protein.Group", "Stripped.Sequence"')
    }
    groups <- DBI::dbGetQuery(con, paste0("SELECT row_group_id, row_group_num_rows, ",
        "stats_min_value, stats_max_value FROM parquet_metadata(?) WHERE path_in_schema = 'Run'"),
        params = list(scope$path))
    counts <- DBI::dbGetQuery(con, 'SELECT "Run", count(*) AS rows FROM qc_input GROUP BY "Run" ORDER BY "Run"')
    stopifnot(nrow(groups) <= 10000L, nrow(counts) == 300L)
    overlap <- vapply(list(samples[1:100], samples[101:200], samples[201:300]), function(selected)
        sum(is.na(groups$stats_min_value) | is.na(groups$stats_max_value) |
            (groups$stats_min_value <= max(selected) & groups$stats_max_value >= min(selected))), integer(1))
    releaseWorkflowTableReadScope(scope)
    stopifnot(manager$.context$store$open_resources == 0L,
        identical(digest::digest(file = original_head, algo = "sha256"), original_hash))
    value <- list(status = "complete", source = input$reference, source_head_sha256 = original_hash,
        scope = "aggregation EXPLAIN plus file/sample metadata; no timed rollup or join attribution",
        row_groups = groups, groups_overlapping_each_batch = overlap,
        sample_rows = counts, plans = plans)
    jsonlite::write_json(value, file.path(root, "rollup-scan-inspection.json"),
        pretty = TRUE, auto_unbox = TRUE, null = "null")
    value
}

#' Frozen native control for the WFX-004 10,000-row filter comparison only.
#' @param manager Existing component owner in the guarded worker.
#' @param reference Frozen input checkpoint.
#' @param parameters Cutoff and inclusion inputs shared with the candidate.
#' @param worker_cap Signature compatibility; the common worker guard still applies.
#' @return The same bounded handoff as the candidate, computed by the native S4 method.
#' @noRd
computeWFXNativeFilterControl <- function(manager, reference, parameters, worker_cap = NULL) {
    before <- readDIAStateMetadata(manager, reference)$data
    scope <- newWorkflowTableReadScope(manager$.context$store,
        before$components$peptide_data, metadata_only = TRUE)
    rows <- scope$metadata$rows
    releaseWorkflowTableReadScope(scope)
    if (rows > 10000L) stop("native comparison control is bounded to the frozen 10,000-row population")
    object <- readDIAExecutionComponents(manager, reference, fields = c("peptide_data", "design_matrix"))
    object <- filterMinNumPeptidesPerSample(object,
        peptides_per_sample_cutoff = parameters$peptides_per_sample_cutoff,
        inclusion_list = parameters$inclusion_list, core_utilisation = NA)
    summary <- object@args$filterMinNumPeptidesPerSample$filter_summary
    # Retain the input design for the common controller's removal bookkeeping.
    checkpoint <- before
    args <- object@args$filterMinNumPeptidesPerSample
    checkpoint$support$filterMinNumPeptidesPerSample <- list(
        support_table = saveDIATableComponent(manager$.context, args$support_table),
        removal_ledger = saveDIATableComponent(manager$.context, args$removal_ledger))
    args$support_table <- args$removal_ledger <- NULL
    checkpoint$slots$args$filterMinNumPeptidesPerSample <- args
    checkpoint$components$peptide_data <- saveDIATableComponent(manager$.context, object@peptide_data)
    checkpoint$dimensions <- list(rows = nrow(object@peptide_data), columns = names(object@peptide_data))
    checkpoint$dependencies <- collectDIANativeReferences(list(reference,
        before$dependencies, before$components, before$support, before$products))
    list(checkpoint = checkpoint, reference = checkpoint$components$peptide_data,
        summary = summary, dependencies = checkpoint$dependencies, counts = NULL,
        execution = list(strategy = "unchanged native S4 filter; fixed 10000-row control"),
        state_name = "sample_filtered", description = "Applied sample quality filter",
        config = parameters, section = "filterMinNumPeptidesPerSample")
}


#' Select an existing replay checkpoint without invoking public S4 hydration.
#' @param manager Owned replay manager; never an original evidence store.
#' @param state_name Exact frozen historical checkpoint name.
#' @return NULL invisibly after reference-only head selection.
#' @noRd
selectDIAReplayCheckpoint <- function(manager, state_name) {
    context <- manager$.context
    snapshot <- context$snapshot
    at <- match(state_name, snapshot$state_history)
    if (is.na(at) || is.null(snapshot$states[[state_name]])) stop("frozen replay state is missing")
    record <- readDIAStateMetadata(manager, snapshot$states[[state_name]])
    if (!isDIAComponentCheckpoint(record$data)) stop("replay input must already be reference-backed")
    snapshot$current_state <- state_name
    snapshot$state_history <- snapshot$state_history[seq_len(at)]
    publishDIAStateHead(context, snapshot)
    invisible(NULL)
}
