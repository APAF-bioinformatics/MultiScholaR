# Fixed PERF-101 diagnostics. Loaded into the retained installed package namespace's
# child environment; this does not replace or trace any scientific function.

#' Run one source-bound, guarded PERF-101 diagnostic.
#' @param library_path,run_dir,repo Explicit retained library and owned paths.
#' @param selection Fixed D0, D1, D0R diagnostic or V verifier selector.
#' @return Diagnostic receipt; never a qualification decision.
#' @noRd
runPerf101Diagnostic <- function(library_path, run_dir, repo, selection) {
    if (identical(selection, "perf101:V")) return(verifyPerf101Outputs(library_path, run_dir, repo))
    attempt <- sub("perf101:", "", selection, fixed = TRUE)
    member <- sub("R$", "", attempt)
    binding <- jsonlite::read_json(file.path(repo,
        "tickets/openai/verification/pilot-source-bindings.json"))
    protocol <- jsonlite::read_json(file.path(repo,
        "tickets/openai/verification/pilot-protocol.json"))
    stopifnot(member %in% c("D0", "D1"),
        identical(normalizePath(dirname(run_dir)), normalizePath(binding$output_root)))
    ledger <- file.path(binding$output_root, paste0(attempt, "-attempt.rds"))
    if (file.exists(ledger)) stop("diagnostic allowance already consumed; no automatic retry")
    if (attempt == "D0R") {
        previous <- readRDS(file.path(binding$output_root, "D0-attempt.rds"))
        failed <- readRDS(file.path(previous$run_dir, "perf101-result.rds"))
        stopifnot(identical(failed$status, "failed"),
            grepl("design_matrix", failed$error, fixed = TRUE))
        correction <- file.path(binding$output_root, "shared-correction.rds")
        if (file.exists(correction)) stop("shared correction pass already consumed")
        saveRDS(list(passes_consumed = 1L, small_rechecks_consumed = 1L,
            remaining_small_rechecks = 1L, owners = c("PERF-101", "PERF-106"),
            reason = "Bind identical design payload across software-specific entity IDs",
            failed_run = previous$run_dir), correction)
    }
    saveRDS(list(member = member, attempt = attempt, run_dir = run_dir, started = Sys.time(),
        deadline_seconds = if (attempt == "D0R") 300 else 600), ledger)
    receipt <- list(status = "running", member = member, protocol = protocol,
        source = binding, phases = list(), qualification = FALSE,
        execution_sources = readRDS(file.path(run_dir, "source-manifest.rds")),
        boundary = "owned physical checkpoint adoption through figures and manager teardown",
        limitations = c("R profiling overhead; diagnostic only", "OS cache uncontrolled",
            "whole service peak includes preparation and independent verification",
            "sampled R stacks do not attribute individual native engine operators"))
    receipt_path <- file.path(run_dir, "perf101-result.rds")
    phase <- function(label, expression) {
        start <- proc.time()
        wall <- Sys.time()
        on.exit({
            receipt$phases[[label]] <<- list(start = wall, end = Sys.time(),
                seconds = unname((proc.time() - start)[["elapsed"]]))
            saveRDS(receipt, receipt_path)
        })
        force(expression)
    }
    on.exit(saveRDS(receipt, receipt_path), add = TRUE)
    tryCatch({
        phase("bindings", {
            frozen <- binding$diagnostic_libraries[[if (member == "D0") "native_D0" else "current_saved_D1"]]
            stopifnot(identical(normalizePath(library_path), normalizePath(frozen$path)))
            check <- function(path, expected) stopifnot(identical(
                digest::digest(file = path, algo = "sha256"), expected))
            check(frozen$rdb_file, frozen$sha256)
            check(binding$source_receipt$path, binding$source_receipt$sha256)
            check(file.path(repo, binding$design_authority$path), binding$design_authority$sha256)
            check(file.path(repo, binding$native_source$path), binding$native_source$sha256)
            check(binding$diagnostic_libraries$authority, binding$diagnostic_libraries$authority_sha256)
            check(binding$peptide_source$path, binding$peptide_source$expected_sha256)
            receipt$session <- utils::sessionInfo()
            receipt$driver_sha256 <- digest::digest(file = file.path(repo,
                "tools/profiling/workflow_perf101.R"), algo = "sha256")
        })
        sys.source(file.path(repo, "tests/testthat/helper-workflow-bulk-rollup.R"),
            envir = environment())
        source_store <- dirname(dirname(dirname(binding$peptide_source$path)))
        source_head <- file.path(source_store, "dia-state-head.rds")
        original_head <- digest::digest(file = source_head, algo = "sha256")
        prepared_path <- file.path(binding$output_root, "screen-input.rds")
        if (member == "D0") {
            phase("completion_failure_observation", {
                results <- lapply(c(0L, 7L), function(status) {
                    out <- processx::run(file.path(R.home("bin"), "Rscript"),
                        c("--vanilla", "-e", paste0("x <- raw(1048576); x[] <- as.raw(1); ",
                            "cat(length(x), '\\n'); quit(status=", status, ")")),
                        error_on_status = FALSE, timeout = 10000)
                    stopifnot(identical(out$status, status), grepl("1048576", out$stdout))
                    out
                })
                receipt$child_observations <- results
                # Both children ran inside this guarded member, including their terminal status.
                # This proves child observation; it is not an OOM/cancellation qualification.
            })
            prepared <- phase("screen_preparation", {
                path <- file.path(run_dir, "screen")
                dir.create(path)
                manager <- copyDIAReplayCheckpoint(source_store, path)
                on.exit(closeDIAStateAdapter(manager), add = TRUE)
                manager$.context$reference_states <- TRUE
                checkpoint <- readDIAStateMetadata(manager)$data
                stopifnot(identical(checkpoint$components$peptide_data$entity_id,
                    binding$peptide_source$entity_id))
                object <- readDIAExecutionComponents(manager, fields = "design_matrix")
                sample <- object@sample_id
                samples <- as.character(object@design_matrix[[sample]])
                authority <- jsonlite::read_json(file.path(repo, binding$design_authority$path))
                predecessor <- readDIAStateMetadata(manager,
                    manager$.context$snapshot$states$protein_peptide_filtered)$data
                design_ref <- predecessor$components$design_matrix
                stopifnot(length(samples) == 300L, !anyNA(samples), !anyDuplicated(samples),
                    identical(design_ref$entity_id, authority$cases$primary$reference$identity$design),
                    identical(design_ref$payload$sha256, checkpoint$components$design_matrix$payload$sha256))
                design <- materializeDIATable(manager$.context, design_ref)
                stopifnot(identical(object@design_matrix, design))
                receipt$design_binding <- list(authority = design_ref,
                    p14 = checkpoint$components$design_matrix, identical_values_types_order = TRUE)
                scope <- newWorkflowTableReadScope(manager$.context$store,
                    checkpoint$components$peptide_data)
                on.exit(releaseWorkflowTableReadScope(scope), add = TRUE)
                con <- scope$connection
                DBI::dbExecute(con, "SET threads=8")
                DBI::dbExecute(con, "SET memory_limit='12GiB'")
                DBI::dbExecute(con, "SET max_temp_directory_size='0B'")
                q <- function(x) as.character(DBI::dbQuoteIdentifier(con, x))
                qs <- function(x) as.character(DBI::dbQuoteString(con, x))
                limits <- rep(c(4000L, 200L, 100L, 4000L), each = 75L)
                rules <- paste0("WHEN ", qs(samples), " THEN ", limits, collapse = " ")
                input <- paste0("read_parquet(", qs(scope$path), ")")
                # The accepted P14 output has one row per protein/stripped-peptide identity.
                # Prove that condition before selecting the first identities by source ordinal.
                counts <- DBI::dbGetQuery(con, paste0("SELECT ", q(sample), " AS sample, count(*) AS n, ",
                    "count(DISTINCT (", q(object@protein_id_column), ",", q(object@peptide_sequence_column),
                    ")) AS identities FROM ", input, " GROUP BY ", q(sample)))
                stopifnot(nrow(counts) == 300L, all(counts$n == 40000),
                    all(counts$n == counts$identities))
                fields <- workflowTableStorageColumns(scope$metadata, scope$metadata$columns)
                key <- if (scope$metadata$key_column == ".rownames" &&
                    isTRUE(scope$metadata$automatic_row_names)) paste0(
                    ', CAST(row_number() OVER (ORDER BY ".multischolar_order") AS VARCHAR)',
                    ' AS ".multischolar_row_key"') else ""
                sql <- paste0("SELECT ", paste(q(fields), collapse = ","), key,
                    ', CAST(row_number() OVER (ORDER BY ".multischolar_order") AS INTEGER)',
                    ' AS ".multischolar_order" FROM (SELECT *, row_number() OVER (PARTITION BY ',
                    q(sample), ' ORDER BY ".multischolar_order") AS perf_rank FROM ', input,
                    ") WHERE perf_rank <= CASE ", q(sample), " ", rules, " ELSE 0 END",
                    ' ORDER BY ".multischolar_order"')
                metadata <- scope$metadata
                metadata$rows <- 622500L
                identity <- newSavedWorkflowScientificIdentity(binding$peptide_source$entity_id,
                    checkpoint$components$design_matrix$entity_id,
                    digest::digest(list(samples, limits), algo = "sha256"),
                    "perf101-screen-v1", manager$.context$software)
                selected <- saveWorkflowTableQuery(scope, sql, metadata, "perf101-screen", identity)
                releaseWorkflowTableReadScope(scope)
                checkpoint$components$peptide_data <- selected
                checkpoint$dimensions <- list(rows = metadata$rows, columns = metadata$columns)
                withDIAStateScope(manager, function() commitDIAReferenceState(manager,
                    "protein_peptide_filtered", checkpoint, list(), "PERF-101 frozen screen"))
                result <- list(store = manager$.context$store$directory, reference = selected,
                    samples = samples, limits = limits, parameters = list(
                        peptides_per_sample_cutoff = 300L, inclusion_list = samples[76:150]),
                    head_sha256 = digest::digest(file = file.path(manager$.context$store$directory,
                        "dia-state-head.rds"), algo = "sha256"))
                closeDIAStateAdapter(manager)
                result
            })
            saveRDS(prepared, prepared_path)
        } else {
            stopifnot(file.exists(prepared_path))
            prepared <- readRDS(prepared_path)
        }
        check(file.path(prepared$store, "dia-state-head.rds"), prepared$head_sha256)
        # Use separate owned stores; keep the prepared input immutable for the paired member.
        owned <- file.path(run_dir, "operation")
        dir.create(owned)
        utils::Rprof(file.path(run_dir, "operation.Rprof"), interval = 0.02, memory.profiling = TRUE)
        utils::Rprofmem(file.path(run_dir, "operation.allocations"), threshold = 1048576)
        on.exit({ utils::Rprof(NULL); utils::Rprofmem(NULL) }, add = TRUE)
        start <- proc.time()[["elapsed"]]
        manager <- phase("adoption", copyDIAReplayCheckpoint(prepared$store, owned))
        on.exit(closeDIAStateAdapter(manager), add = TRUE)
        manager$.context$reference_states <- TRUE
        value <- phase("load_validate_compute_write", {
            if (member == "D0") computePerf101Native(manager, getDIAStateReference(manager),
                prepared$parameters) else computeDIAReferenceSampleFilter(manager,
                    getDIAStateReference(manager), prepared$parameters)
        })
        receipt$inline_result_bytes <- as.numeric(object.size(value))
        receipt$summary <- value$summary
        phase("summary_checkpoint_publication", {
            workflow <- new.env(parent = emptyenv())
            workflow$state_manager <- manager
            workflow$qc_params <- list(peptide_qc = list())
            product <- withDIAStateScope(manager, function() finishDIAReferencePeptideStep(workflow,
                value$checkpoint, value$state_name, value$description, value$config, value$section,
                value$counts, observed_inclusion_absence = value$summary$inclusion_entries_absent_from_data))
            owner <- list(manager = manager)
            publishWorkflowPeptideResult(owner, list(request_dir = run_dir),
                list(state_reference = product$stateReference, operation = "minimum_peptides_per_sample"))
            assertDIAExportReady(manager)
            receipt$state <- product$stateReference
        })
        phase("teardown", {
            stopifnot(manager$.context$store$open_resources == 0L)
            closeDIAStateAdapter(manager)
        })
        receipt$operation_seconds <- proc.time()[["elapsed"]] - start
        utils::Rprof(NULL)
        utils::Rprofmem(NULL)
        receipt$profile <- utils::summaryRprof(file.path(run_dir, "operation.Rprof"), memory = "both")
        phase("independent_verification", {
            manager <- newDIAStateAdapter(newWorkflowEntityStore(owned,
                prepared$reference$owner$project_id, prepared$reference$owner$workflow_id), reopen = TRUE)
            scope <- newWorkflowTableReadScope(manager$.context$store, value$reference, reader = "fastduckdb")
            on.exit(releaseWorkflowTableReadScope(scope), add = TRUE)
            control <- if (member == "D1") bindPerf101Control(binding$output_root) else NULL
            comparison <- NULL
            if (!is.null(control)) {
                comparison <- openPerf101Comparison(control)
                on.exit({
                    releaseWorkflowTableReadScope(comparison$scope)
                    closeWorkflowEntityStore(comparison$store)
                }, add = TRUE)
                receipt$native_control <- control[c("receipt_path", "receipt_sha256")]
                stopifnot(comparison$scope$metadata$rows == scope$metadata$rows,
                    identical(comparison$scope$metadata$columns, scope$metadata$columns))
            }
            runs <- character()
            read_page <- function(offset, limit) {
                page <- readSavedWorkflowTable(manager$.context$store, value$reference,
                    scope$metadata$columns, offset = offset, limit = limit, scope = scope)
                runs <<- union(runs, as.character(page$data$Run))
                page
            }
            other_page <- if (is.null(comparison)) NULL else function(offset, limit)
                readSavedWorkflowTable(comparison$store, comparison$scope$reference,
                    comparison$scope$metadata$columns, offset = offset, limit = limit,
                    scope = comparison$scope)
            hashes <- comparePerf101Pages(read_page, other_page, scope$metadata$rows)
            rows <- scope$metadata$rows
            releaseWorkflowTableReadScope(scope)
            stopifnot(rows == 615000L, setequal(runs, prepared$samples[-(151:225)]))
            receipt$output_hashes <- hashes
            receipt$output_rows <- rows
            receipt$output_samples <- length(runs)
            receipt$support <- lapply(value$checkpoint$support$filterMinNumPeptidesPerSample,
                function(ref) materializeDIATable(manager$.context, ref))
            figures <- list.files(file.path(run_dir, "required-output"),
                pattern = "[.](png|pdf)$", recursive = TRUE, full.names = TRUE)
            receipt$figures <- figures
            stopifnot(length(figures) >= 12L, all(file.size(figures) > 0))
            if (member == "D1") {
                stopifnot(identical(receipt$summary, control$result$summary),
                    identical(receipt$support, control$result$support))
                releaseWorkflowTableReadScope(comparison$scope)
                closeWorkflowEntityStore(comparison$store)
            }
            closeDIAStateAdapter(manager)
            check(source_head, original_head)
            check(binding$peptide_source$path, binding$peptide_source$expected_sha256)
            check(file.path(prepared$store, "dia-state-head.rds"), prepared$head_sha256)
            receipt$science <- TRUE
            receipt$original_unchanged <- TRUE
        })
        receipt$status <- "complete"
        saveRDS(receipt, receipt_path)
        receipt
    }, error = function(error) {
        receipt$status <<- "failed"
        receipt$error <<- conditionMessage(error)
        saveRDS(receipt, receipt_path)
        stop(error)
    })
}

#' Adapt the existing native control's persistence to the frozen PERF-101 population.
#' The scientific S4 method is unchanged and called exactly once.
#' @param manager,reference,parameters Existing execution inputs.
#' @return Existing compute-only checkpoint envelope.
#' @noRd
computePerf101Native <- function(manager, reference, parameters) {
    checkpoint <- readDIAStateMetadata(manager, reference)$data
    object <- readDIAExecutionComponents(manager, reference, fields = c("peptide_data", "design_matrix"))
    stopifnot(nrow(object@peptide_data) == 622500L)
    object <- filterMinNumPeptidesPerSample(object,
        peptides_per_sample_cutoff = parameters$peptides_per_sample_cutoff,
        inclusion_list = parameters$inclusion_list, core_utilisation = NA)
    args <- object@args$filterMinNumPeptidesPerSample
    summary <- args$filter_summary
    checkpoint$support$filterMinNumPeptidesPerSample <- list(
        support_table = saveDIATableComponent(manager$.context, args$support_table),
        removal_ledger = saveDIATableComponent(manager$.context, args$removal_ledger))
    args$support_table <- args$removal_ledger <- NULL
    checkpoint$slots$args$filterMinNumPeptidesPerSample <- args
    checkpoint$components$peptide_data <- saveDIATableComponent(manager$.context, object@peptide_data)
    checkpoint$dimensions <- list(rows = nrow(object@peptide_data), columns = names(object@peptide_data))
    checkpoint$dependencies <- collectDIANativeReferences(list(reference, checkpoint$dependencies,
        checkpoint$components, checkpoint$support, checkpoint$products))
    list(checkpoint = checkpoint, reference = checkpoint$components$peptide_data,
        summary = summary, dependencies = checkpoint$dependencies, counts = NULL,
        execution = list(strategy = "unchanged native S4 filter; single process"),
        state_name = "sample_filtered", description = "Applied sample quality filter",
        config = parameters, section = "filterMinNumPeptidesPerSample")
}

#' Compare already-produced diagnostic outputs using fixed logical row pages.
#' @param library_path,run_dir,repo Explicit retained library and owned paths.
#' @return A separate verification receipt; original failed receipts stay unchanged.
#' @noRd
verifyPerf101Outputs <- function(library_path, run_dir, repo) {
    binding <- jsonlite::read_json(file.path(repo,
        "tickets/openai/verification/pilot-source-bindings.json"))
    root <- binding$output_root
    stopifnot(identical(normalizePath(dirname(run_dir)), normalizePath(root)))
    ledger <- file.path(root, "V-attempt.rds")
    if (file.exists(ledger)) stop("verifier-only allowance already consumed")
    saveRDS(list(run_dir = run_dir, started = Sys.time(), seconds = 300L,
        authority = "user continue after explicit verifier-only correction request",
        scientific_reexecution = FALSE, additional_correction_pass = 1L), ledger)
    receipt <- list(status = "running", filtering_rerun = FALSE, qualified = FALSE,
        execution_sources = readRDS(file.path(run_dir, "source-manifest.rds")),
        page_limit = 10000L, completed_pages = 0L)
    path <- file.path(run_dir, "perf101-verification.rds")
    started <- proc.time()[["elapsed"]]
    on.exit({
        receipt$elapsed_seconds <- proc.time()[["elapsed"]] - started
        saveRDS(receipt, path)
    }, add = TRUE)
    tryCatch({
        frozen <- binding$diagnostic_libraries$current_saved_D1
        stopifnot(identical(normalizePath(library_path), normalizePath(frozen$path)),
            identical(digest::digest(file = frozen$rdb_file, algo = "sha256"), frozen$sha256))
        control <- bindPerf101Control(root)
        member_dir <- readRDS(file.path(root, "D1-attempt.rds"))$run_dir
        result_path <- file.path(member_dir, "perf101-result.rds")
        candidate <- list(directory = member_dir, result = readRDS(result_path),
            receipt_path = result_path, receipt_sha256 = digest::digest(file = result_path, algo = "sha256"))
        runs <- list(control, candidate)
        # A failed post-compute verification is eligible only if all operation phases
        # and the durable published state exist. There is no required error wording.
        stopifnot(candidate$result$status %in% c("complete", "failed"),
            !is.null(candidate$result$state), !is.null(candidate$result$phases$teardown),
            identical(control$result$summary, candidate$result$summary),
            identical(control$result$support, candidate$result$support))
        receipt$source_receipts <- lapply(runs, function(x) x[c("receipt_path", "receipt_sha256")])
        stores <- list()
        scopes <- list()
        on.exit({
            for (scope in scopes) releaseWorkflowTableReadScope(scope)
            for (store in stores) closeWorkflowEntityStore(store)
        }, add = TRUE)
        checkpoints <- list()
        for (i in seq_along(runs)) {
            project <- file.path(runs[[i]]$directory, "operation")
            owner <- readSavedStoreIdentity(file.path(project, ".workflow-entities/store.json"))
            stores[[i]] <- newWorkflowEntityStore(project, owner$project_id, owner$workflow_id)
            state <- runs[[i]]$result$state
            validateSavedWorkflowReference(stores[[i]], state, payload = TRUE)
            checkpoints[[i]] <- readRDS(file.path(stores[[i]]$directory,
                state$payload$relative_path))$data
            scopes[[i]] <- newWorkflowTableReadScope(stores[[i]],
                checkpoints[[i]]$components$peptide_data, reader = "fastduckdb")
            stopifnot(scopes[[i]]$metadata$rows == 615000L)
            figures <- runs[[i]]$result$figures
            stopifnot(length(figures) == 12L, all(file.exists(figures)), all(file.size(figures) > 0))
        }
        stopifnot(identical(scopes[[1L]]$metadata$columns, scopes[[2L]]$metadata$columns))
        fetch <- function(i) function(offset, limit) readSavedWorkflowTable(stores[[i]],
            scopes[[i]]$reference, scopes[[i]]$metadata$columns, offset = offset,
            limit = limit, scope = scopes[[i]], reader = "fastduckdb")
        receipt$block_hashes <- comparePerf101Pages(fetch(1L), fetch(2L), 615000L)
        receipt$completed_pages <- length(receipt$block_hashes)
        designs <- lapply(seq_along(stores), function(i) {
            ref <- checkpoints[[i]]$components$design_matrix
            readSavedWorkflowTable(stores[[i]], ref,
                columns = readSavedWorkflowTableMetadata(file.path(stores[[i]]$directory,
                    ref$payload$relative_path))$columns, limit = 1000L, reader = "fastduckdb")$data
        })
        stopifnot(identical(designs[[1L]], designs[[2L]]), nrow(designs[[1L]]) == 225L)
        for (i in seq_along(scopes)) releaseWorkflowTableReadScope(scopes[[i]])
        stopifnot(all(vapply(stores, function(store) store$open_resources == 0L, logical(1))))
        for (item in runs) stopifnot(identical(digest::digest(file = item$receipt_path,
            algo = "sha256"), item$receipt_sha256))
        stopifnot(identical(digest::digest(file = binding$peptide_source$path, algo = "sha256"),
            binding$peptide_source$expected_sha256))
        receipt$rows <- 615000L
        receipt$samples <- 225L
        receipt$exact_values_attributes_order <- TRUE
        receipt$design_summary_support_figures <- TRUE
        receipt$resources_released <- TRUE
        receipt$original_receipts_unchanged <- TRUE
        receipt$status <- "complete"
        receipt
    }, error = function(error) {
        receipt$status <<- "failed"
        receipt$error <<- conditionMessage(error)
        stop(error)
    })
}

#' Compare fixed logical pages, independent of physical Parquet batch boundaries.
#' @param left,right Functions taking offset/limit and returning existing codec pages.
#' @param rows Exact logical row count; right may be NULL for the native oracle.
#' @return Ordered hashes including data attributes, row names and source order.
#' @noRd
comparePerf101Pages <- function(left, right = NULL, rows) {
    stopifnot(is.function(left), is.null(right) || is.function(right),
        length(rows) == 1L, is.finite(rows), rows >= 0, rows == floor(rows))
    hashes <- character(ceiling(rows / 10000))
    for (i in seq_along(hashes)) {
        offset <- (i - 1L) * 10000L
        limit <- min(10000L, rows - offset)
        page <- left(offset, limit)
        stopifnot(is.data.frame(page$data), nrow(page$data) == limit,
            length(page$source_order) == limit, length(page$source_row_names) == limit)
        fields <- c("data", "source_order", "source_row_names")
        if (!is.null(right)) {
            other <- right(offset, limit)
            if (!identical(page[fields], other[fields]))
                stop("logical page differs at row ", offset + 1L)
        }
        hashes[[i]] <- digest::digest(page[fields], algo = "sha256")
    }
    hashes
}

#' Resolve and bind the one completed native control; never assume a retry.
#' @param root Frozen diagnostic evidence root.
#' @return Explicit receipt path, hash, directory and verified completed result.
#' @noRd
bindPerf101Control <- function(root) {
    controls <- list()
    for (member in c("D0", "D0R")) {
        attempt <- file.path(root, paste0(member, "-attempt.rds"))
        if (!file.exists(attempt)) next
        directory <- readRDS(attempt)$run_dir
        path <- file.path(directory, "perf101-result.rds")
        if (!file.exists(path)) next
        result <- readRDS(path)
        if (!identical(result$status, "complete")) next
        stopifnot(identical(result$member, "D0"))
        controls[[length(controls) + 1L]] <- list(directory = directory,
            receipt_path = path, receipt_sha256 = digest::digest(file = path, algo = "sha256"),
            result = result)
    }
    if (length(controls) != 1L) stop("expected exactly one completed native control")
    controls[[1L]]
}

#' Open a bound native result with the existing codec for page comparison.
#' @param control Immutable result/receipt binding from bindPerf101Control.
#' @return Store and scoped table owner, closed by the caller.
#' @noRd
openPerf101Comparison <- function(control) {
    root <- file.path(control$directory, "operation")
    owner <- readSavedStoreIdentity(file.path(root, ".workflow-entities/store.json"))
    store <- newWorkflowEntityStore(root, owner$project_id, owner$workflow_id)
    ready <- FALSE
    on.exit(if (!ready) closeWorkflowEntityStore(store))
    state <- control$result$state
    validateSavedWorkflowReference(store, state, payload = TRUE)
    checkpoint <- readRDS(file.path(store$directory, state$payload$relative_path))$data
    scope <- newWorkflowTableReadScope(store, checkpoint$components$peptide_data,
        reader = "fastduckdb")
    ready <- TRUE
    list(store = store, scope = scope)
}
