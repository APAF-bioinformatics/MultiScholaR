# Fixed PERF-106 qualification through the existing scoped codecs and runner.

#' Resolve scientific verification only after its owned service finishes cleanly.
#' @param owned Operation member directory containing the verification launch log.
#' @param output_root Existing service output root.
#' @param exit_status Actual verifier launcher exit status, never inferred from a receipt.
#' @return Accepted bound value/terminal identity, or an explicit rejection reason.
#' @noRd
readPerf106Verification <- function(owned, output_root, exit_status) {
    reject <- function(reason) list(accepted = FALSE, reason = reason, value = NULL)
    if (length(exit_status) != 1L || !is.numeric(exit_status) || is.na(exit_status) || exit_status != 0)
        return(reject("verifier launcher did not exit successfully"))
    tryCatch({
        log <- file.path(owned, "verification-launch.log")
        path <- file.path(owned, "verification.rds")
        operation <- file.path(owned, "operation-receipt.rds")
        if (any(!file.exists(c(log, path, operation))) || file.size(log) > 1048576 ||
            file.size(path) > 8388608) return(reject("verification evidence missing or oversized"))
        launches <- grep('^\\{"unit":', readLines(log, warn = FALSE), value = TRUE)
        if (length(launches) != 1L) return(reject("verifier launch identity is ambiguous"))
        launch <- jsonlite::fromJSON(launches)
        if (!is.character(launch$unit) || length(launch$unit) != 1L ||
            !grepl("^multischolar-dia003-browser-smoke[.][A-Za-z0-9]{6}[.]service$", launch$unit) ||
            !is.character(launch$invocation_id) || length(launch$invocation_id) != 1L ||
            !grepl("^[a-f0-9]{32}$", launch$invocation_id)) return(reject("invalid verifier launch identity"))
        directory <- file.path(normalizePath(output_root),
            sub("[.]service$", "", sub("^multischolar-dia003-", "", launch$unit)))
        terminal_path <- file.path(directory, "terminal-interval.json")
        controller_path <- file.path(directory, "controller-status.txt")
        if (!file.exists(terminal_path) || !file.exists(controller_path) ||
            file.exists(file.path(directory, "guard-stop.txt"))) return(reject("verifier terminal/guard evidence failed"))
        terminal <- jsonlite::read_json(terminal_path)
        controller <- readLines(controller_path, warn = FALSE)
        expected_controller <- c(paste0("unit=", sub("[.]service$", "", launch$unit)),
            "interval_exit=0", "stop_exit=0")
        value <- readRDS(path)
        if (!identical(terminal$status, "service_tree_interval_only") ||
            !isTRUE(terminal$service_succeeded) || !identical(terminal$service_result, "success") ||
            !identical(terminal$unit, launch$unit) || !identical(terminal$invocation_id, launch$invocation_id) ||
            !identical(sort(controller), sort(expected_controller)) ||
            !identical(value$status, "complete") ||
            !identical(normalizePath(value$run_dir), normalizePath(directory)) ||
            !identical(value$operation_sha256, digest::digest(file = operation, algo = "sha256")))
            return(reject("verifier completion is unsuccessful or differs from the operation"))
        list(accepted = TRUE, reason = NULL, value = value, unit = launch$unit,
            invocation_id = launch$invocation_id,
            terminal_sha256 = digest::digest(file = terminal_path, algo = "sha256"),
            verification_sha256 = digest::digest(file = path, algo = "sha256"))
    }, error = function(error) reject(paste("invalid verifier evidence:", conditionMessage(error))))
}

#' Build an untimed worker-entry profiler without wrapping scientific arguments.
#' @param entry Fixed installed entry, or the instrumentation check's fixed fixture.
#' @return Rscript expression writing profiles beside its owned request.
#' @noRd
perf106ProfileExpression <- function(entry = "MultiScholaR:::workflowWorkerMain(commandArgs(TRUE))") {
    paste0('local({ d <- dirname(commandArgs(TRUE)[[1L]]); ',
        'clock <- function() as.numeric(strsplit(readLines("/proc/uptime", n=1L), " ")[[1L]][[1L]]); ',
        'saveRDS(list(start=clock(), pid=Sys.getpid()), file.path(d,"profile-start.rds")); ',
        'utils::Rprof(file.path(d,"worker-stack.out"), interval=0.01, memory.profiling=TRUE, ',
        'gc.profiling=TRUE, event="elapsed"); ',
        'utils::Rprofmem(file.path(d,"worker-allocations.out"), threshold=1048576); ',
        'failed <- FALSE; tryCatch({ ', entry, ' }, error=function(e) { failed <<- TRUE; stop(e) }, ',
        'finally={ utils::Rprofmem(NULL); utils::Rprof(NULL); ',
        'saveRDS(list(finished=clock(), error=failed), file.path(d,"profile-finished.rds")) }) })')
}

#' Replace only the fixed worker command in a private diagnostic launcher copy.
#' @param launch Existing installed service launcher.
#' @return Same launcher with profiling at process entry; guards/argv remain intact.
#' @noRd
perf106ProfileLauncher <- function(launch, expression = perf106ProfileExpression()) {
    original <- "MultiScholaR:::workflowWorkerMain(commandArgs(TRUE))"
    replacement <- expression
    count <- 0L
    rewrite <- function(value) {
        if (identical(value, original)) { count <<- count + 1L; return(replacement) }
        if (is.call(value)) return(as.call(lapply(as.list(value), rewrite)))
        value
    }
    body(launch) <- rewrite(body(launch))
    stopifnot(count == 1L)
    launch
}

#' Prove profiler completion/failure flushing and shared-clock observations.
#' @param run_dir,repo Guarded instrumentation output and frozen source.
#' @return Fixture-only evidence, never scientific qualification.
#' @noRd
checkPerf106ProfileInstrumentation <- function(run_dir, repo) {
    binding <- jsonlite::read_json(file.path(repo,
        "tickets/openai/verification/pilot-source-bindings.json"))
    ledger <- file.path(binding$output_root, "repair-v2/profile-check.rds")
    stopifnot(!file.exists(ledger))
    result <- list(status = "starting", run_dir = run_dir)
    saveRDS(result, ledger)
    tryCatch({
        # Only the launch expression changes. The installed worker and its
        # scientific functions are neither copied nor rebound in its process.
        launch <- perf106ProfileLauncher(workflowLaunchService)
        stopifnot(identical(formals(launch), formals(workflowLaunchService)),
            identical(environment(launch), environment(workflowLaunchService)))
        expression <- perf106ProfileExpression(paste(
            'x <- rep.int(as.raw(1L), 33554432L); Sys.sleep(0.2);',
            'if (isTRUE(readRDS(commandArgs(TRUE)[[1L]])$fail)) stop("owned profile failure");',
            'stopifnot(length(x)==33554432L)'))
        result$cases <- list()
        for (failed in c(FALSE, TRUE)) {
            owned <- file.path(run_dir, if (failed) "failure" else "completion")
            stopifnot(dir.create(owned))
            request <- file.path(owned, "request.rds")
            saveRDS(list(fail = failed), request)
            before <- as.numeric(strsplit(readLines("/proc/uptime", n = 1L), " ")[[1L]][[1L]])
            status <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
                c("--vanilla", "-e", shQuote(expression), shQuote(request)),
                stdout = file.path(owned, "stdout.log"), stderr = file.path(owned, "stderr.log")))
            after <- as.numeric(strsplit(readLines("/proc/uptime", n = 1L), " ")[[1L]][[1L]])
            start <- readRDS(file.path(owned, "profile-start.rds"))
            finish <- readRDS(file.path(owned, "profile-finished.rds"))
            stack <- readLines(file.path(owned, "worker-stack.out"))
            allocations <- readLines(file.path(owned, "worker-allocations.out"))
            sizes <- suppressWarnings(as.numeric(sub(" .*", "", allocations)))
            stopifnot(identical(status == 0L, !failed), identical(finish$error, failed),
                start$start >= before, finish$finished <= after,
                finish$finished > start$start, length(stack) > 2L,
                any(sizes >= 33554432, na.rm = TRUE))
            result$cases[[basename(owned)]] <- list(exit_status = status,
                start = start, finish = finish, largest_R_allocation = max(sizes, na.rm = TRUE))
        }
        result$status <- "complete"
    }, error = function(error) { result$status <<- "failed"; result$error <<- conditionMessage(error) })
    saveRDS(result, ledger)
    saveRDS(result, file.path(run_dir, "profile-check.rds"))
    if (result$status != "complete") stop(result$error)
    result
}

#' Verify the actual P14 executor against the immutable accepted rollup oracle.
#' @param run_dir,repo Guarded verification output and frozen source snapshot.
#' @param slot Fixed P14 verification index, twenty-one through twenty-four.
#' @return Independent science, figures and execution identity receipt.
#' @noRd
verifyPerf106RepairP14 <- function(run_dir, repo, slot, stage = "p14") {
    stopifnot(stage %in% c("p14", "p14-diagnostic"))
    stopifnot(slot %in% 21:24)
    source <- jsonlite::read_json(file.path(repo,
        "tickets/openai/verification/pilot-source-bindings.json"))
    root <- file.path(source$output_root, "repair-v2")
    binding <- readRDS(file.path(repo,
        "docs/qa/workflow-execution/perf-106/p14-regression-bindings.rds"))
    i <- slot - 20L
    role <- c("control", "candidate", "candidate", "control")[[i]]
    rep <- c(1L, 1L, 2L, 2L)[[i]]
    owned <- file.path(root, stage, paste(role, rep, sep = "-"))
    ledger <- file.path(owned, "verification.rds")
    stopifnot(!file.exists(ledger))
    receipt <- list(status = "starting", run_dir = run_dir, role = role, rep = rep)
    saveRDS(receipt, ledger)
    tryCatch({
        operation_path <- file.path(owned, "operation-receipt.rds")
        receipt$operation_sha256 <- digest::digest(file = operation_path, algo = "sha256")
        operation <- readRDS(operation_path)
        member <- readRDS(file.path(owned, "member-result.rds"))
        builds <- readRDS(file.path(root, stage, "stage.rds"))$p14_builds
        stopifnot(identical(operation$status, "complete"), identical(member$status, "complete"),
            identical(operation$config$build, builds[[role]]),
            identical(member$project_root, file.path(owned, "project")),
            identical(digest::digest(file = binding$accepted_receipt, algo = "sha256"),
                binding$accepted_receipt_sha256))
        project <- member$project_root
        owner <- readSavedStoreIdentity(file.path(project, ".workflow-entities/store.json"))
        manager <- newDIAStateAdapter(newWorkflowEntityStore(project,
            owner$project_id, owner$workflow_id), reopen = TRUE)
        on.exit(closeDIAStateAdapter(manager), add = TRUE)
        checkpoint <- readDIAStateMetadata(manager)$data
        reference <- checkpoint$components$peptide_data
        sys.source(file.path(repo, "tests/testthat/helper-workflow-bulk-rollup.R"),
            envir = environment())
        hashes <- fingerprintDIAWriterOutput(project, reference)
        encode <- function(x) jsonlite::toJSON(x, auto_unbox = TRUE, null = "null", digits = NA)
        request <- file.path(project, ".workflow-entities/.workflow-export",
            operation$execution$request_id)
        worker <- readRDS(file.path(request, "worker-result.rds"))
        execution <- worker$result$execution
        figures <- list.files(file.path(request, "required-output"),
            "4_precursor_rollup_.*[.](png|pdf)$", recursive = TRUE, full.names = TRUE)
        figures <- figures[!grepl("/figures/", figures)]
        stopifnot(reference$metadata$shape[[1L]] == 12000000L, length(hashes) == 1200L,
            identical(unname(hashes), unname(binding$accepted_output_hashes)),
            identical(encode(checkpoint$slots$args$rollUpPrecursorToPeptide$rollup_summary),
                encode(binding$accepted_summary)),
            identical(execution$exact_sum_sha256, binding$accepted_exact_sum_sha256),
            identical(execution$profile$name, "P14"), identical(execution$engine_threads, 8L),
            identical(execution$writer_backend, "fast_r"), identical(execution$writer_threads, 8L),
            !is.null(execution$admission$worker_cap),
            length(figures) == 12L, all(file.size(figures) > 0),
            is.null(assertDIAExportReady(manager)),
            identical(digest::digest(file = file.path(binding$source_store, "dia-state-head.rds"),
                algo = "sha256"), binding$source_head_sha256),
            identical(digest::digest(file = operation_path, algo = "sha256"), receipt$operation_sha256))
        receipt$status <- "complete"
        receipt$pages <- hashes
        receipt$summary_exact <- TRUE
        receipt$figures <- stats::setNames(file.size(figures), basename(figures))
        receipt$execution <- execution
    }, error = function(error) {
        receipt$status <<- "failed"
        receipt$error <<- conditionMessage(error)
    })
    saveRDS(receipt, ledger)
    saveRDS(receipt, file.path(run_dir, "perf106-p14-verification.rds"))
    if (receipt$status != "complete") stop(receipt$error)
    receipt
}

#' Verify one frozen screen output outside its measured application slice.
#' @param library_path,run_dir,repo Guarded verifier build, output and source snapshot.
#' @param slot Fixed member index: screen one through twelve; confirmation thirteen through twenty.
#' @return Immutable verification receipt bound to the operation receipt.
#' @noRd
verifyPerf106RepairScreen <- function(library_path, run_dir, repo, slot) {
    stopifnot(slot %in% seq_len(20L))
    binding <- jsonlite::read_json(file.path(repo,
        "tickets/openai/verification/pilot-source-bindings.json"))
    root <- file.path(binding$output_root, "repair-v2")
    input <- readRDS(file.path(root, "qualification-inputs.rds"))
    full <- slot > 12L
    stage <- if (full) "confirmation" else "screen"
    local_slot <- if (full) slot - 12L else slot
    width <- if (full) 4L else 6L
    population <- c("primary", "wide")[[ceiling(local_slot / width)]]
    index <- (local_slot - 1L) %% width + 1L
    if (full) {
        decision <- readRDS(file.path(root, "screen/decision.rds"))
        stopifnot(identical(decision$status, "screen_finalist"))
        arm <- c("native", decision$finalist, decision$finalist, "native")[[index]]
        rep <- c(1L, 1L, 2L, 2L)[[index]]
    } else {
        arm <- c("native", "physical", "selection", "selection", "physical", "native")[[index]]
        rep <- c(1L, 1L, 1L, 2L, 2L, 2L)[[index]]
    }
    id <- paste(population, arm, rep, sep = "-")
    owned <- file.path(root, stage, id)
    ledger <- file.path(owned, "verification.rds")
    stopifnot(!file.exists(ledger))
    receipt <- list(status = "starting", run_dir = run_dir, member_id = id)
    saveRDS(receipt, ledger)
    tryCatch({
        operation_path <- file.path(owned, "operation-receipt.rds")
        receipt$operation_sha256 <- digest::digest(file = operation_path, algo = "sha256")
        operation <- readRDS(operation_path)
        member <- readRDS(file.path(owned, "member-result.rds"))
        expected <- input[[if (full) "full" else "screen"]][[population]]
        stopifnot(identical(operation$status, "complete"),
            identical(member$status, "complete"),
            identical(operation$config$build, input$build),
            identical(workflowExecutionBuildIdentity(), input$build),
            identical(operation$member$parameters, expected$parameters),
            identical(operation$member$input_entity_id, expected$reference$entity_id),
            identical(member$project_root, file.path(owned, "project")))
        project <- member$project_root
        owner <- readSavedStoreIdentity(file.path(project, ".workflow-entities/store.json"))
        manager <- newDIAStateAdapter(newWorkflowEntityStore(project,
            owner$project_id, owner$workflow_id), reopen = TRUE)
        on.exit(closeDIAStateAdapter(manager), add = TRUE)
        reference <- getDIAStateReference(manager)
        receipt$verification <- verifyPerf106Result(project, reference)
        worker_path <- file.path(project, ".workflow-entities", ".workflow-export",
            operation$execution$request_id, "worker-result.rds")
        worker <- readRDS(worker_path)
        receipt$science_execution <- worker$result$execution
        stopifnot(identical(receipt$science_execution$bulk_query$settings,
            list(threads = 8L, engine_limit = "12GiB", spill = FALSE)))
        stopifnot(receipt$verification$rows == if (full) 6150000L else 615000L,
            receipt$verification$samples == 225L,
            isTRUE(receipt$verification$required_outputs_complete),
            identical(digest::digest(file = operation_path, algo = "sha256"), receipt$operation_sha256),
            identical(digest::digest(file = file.path(expected$store, "dia-state-head.rds"),
                algo = "sha256"), expected$head_sha256),
            identical(digest::digest(file = file.path(expected$store,
                expected$reference$payload$relative_path), algo = "sha256"), expected$payload_sha256))
        receipt$status <- "complete"
    }, error = function(error) {
        receipt$status <<- "failed"
        receipt$error <<- conditionMessage(error)
    })
    saveRDS(receipt, ledger)
    saveRDS(receipt, file.path(run_dir, "perf106-screen-verification.rds"))
    if (receipt$status != "complete") stop(receipt$error)
    receipt
}

#' Freeze retained populations and prepare only the missing wide screen input.
#' @param library_path,run_dir,repo Installed package, guarded output and source snapshot.
#' @return Input/build bindings and exact bounded-page preparation evidence.
#' @noRd
preparePerf106RepairInputs <- function(library_path, run_dir, repo) {
    binding <- jsonlite::read_json(file.path(repo,
        "tickets/openai/verification/pilot-source-bindings.json"))
    root <- file.path(binding$output_root, "repair-v2")
    ledger <- file.path(root, "input-preparation-attempt.rds")
    output <- file.path(root, "qualification-inputs.rds")
    stopifnot(!file.exists(ledger), !file.exists(output))
    receipt <- list(status = "starting", started = Sys.time(), run_dir = run_dir)
    flush <- function() {
        saveRDS(receipt, ledger)
        saveRDS(receipt, file.path(run_dir, "perf106-repair-inputs.rds"))
    }
    flush()
    tryCatch({
        sys.source(file.path(repo, "tests/testthat/helper-workflow-bulk-rollup.R"),
            envir = environment())
        validate <- function(input, rows) {
            stopifnot(identical(digest::digest(file = file.path(input$store,
                "dia-state-head.rds"), algo = "sha256"), input$head_sha256))
            owner <- readSavedStoreIdentity(file.path(input$store, "store.json"))
            manager <- newDIAStateAdapter(newWorkflowEntityStore(dirname(input$store),
                owner$project_id, owner$workflow_id), reopen = TRUE)
            on.exit(closeDIAStateAdapter(manager))
            checkpoint <- readDIAStateMetadata(manager)$data
            stopifnot(identical(checkpoint$components$peptide_data, input$reference))
            validateSavedWorkflowReference(manager$.context$store, input$reference, payload = TRUE)
            payload <- file.path(input$store, input$reference$payload$relative_path)
            hash <- digest::digest(file = payload, algo = "sha256")
            if (!is.null(input$payload_sha256)) stopifnot(identical(hash, input$payload_sha256))
            stopifnot(input$reference$metadata$shape[[1L]] == rows,
                length(input$samples) == 300L, !anyDuplicated(input$samples))
            design <- materializeDIATable(manager$.context, checkpoint$components$design_matrix)
            stopifnot(identical(as.character(design[[checkpoint$slots$sample_id]]), input$samples))
            input$payload_sha256 <- hash
            input$payload_bytes <- as.numeric(file.size(payload))
            input$rows <- as.integer(rows)
            input
        }
        full <- readRDS(file.path(repo,
            "docs/qa/workflow-execution/perf-106/confirmation-input-bindings.rds"))
        full <- lapply(full, validate, rows = 6225000L)
        primary <- validate(readRDS(file.path(binding$output_root, "screen-input.rds")), 622500L)
        primary$annotation_suffix_bytes <- 0L
        owned <- file.path(run_dir, "wide-screen")
        stopifnot(dir.create(owned))
        manager <- copyDIAReplayCheckpoint(primary$store, owned)
        on.exit(closeDIAStateAdapter(manager), add = TRUE)
        manager$.context$reference_states <- TRUE
        checkpoint <- readDIAStateMetadata(manager)$data
        scope <- newWorkflowTableReadScope(manager$.context$store,
            checkpoint$components$peptide_data, reader = "fastduckdb")
        scope$read_threads <- 8L
        on.exit(releaseWorkflowTableReadScope(scope), add = TRUE, after = FALSE)
        metadata <- scope$metadata
        stopifnot(identical(metadata$column_metadata$Protein.Ids$typeof, "character"))
        suffix <- strrep("x", 2048L)
        identity <- newSavedWorkflowScientificIdentity(primary$reference$entity_id,
            checkpoint$components$design_matrix$entity_id,
            digest::digest(list(primary$parameters, suffix), algo = "sha256"),
            "perf106-repair-wide-screen", manager$.context$software)
        wide <- saveWorkflowTable(scope$store, NULL, "wide-screen", identity,
            metadata$key_column, metadata = metadata, write_fn = function(path, restoration) {
                writer <- openWorkflowTableWriter(path, restoration, backend = "fast_r")
                on.exit(abortWorkflowTableWriter(writer), add = TRUE)
                appendWorkflowTableWriter(writer, readSavedWorkflowTablePrototype(scope))
                streamSavedWorkflowTable(scope, metadata$columns, function(data, offset, row_names) {
                    data$Protein.Ids <- paste0(data$Protein.Ids, suffix)
                    appendWorkflowTableWriter(writer, data)
                })
                finishWorkflowTableWriter(writer)
            })
        wide_scope <- newWorkflowTableReadScope(scope$store, wide, reader = "fastduckdb")
        wide_scope$read_threads <- 8L
        on.exit(releaseWorkflowTableReadScope(wide_scope), add = TRUE, after = FALSE)
        page_hashes <- character(ceiling(metadata$rows / 10000L))
        for (i in seq_along(page_hashes)) {
            offset <- (i - 1L) * 10000L
            limit <- min(10000L, metadata$rows - offset)
            left <- readSavedWorkflowTable(scope$store, scope$reference, metadata$columns,
                offset = offset, limit = limit, scope = scope)$data
            right <- readSavedWorkflowTable(scope$store, wide, metadata$columns,
                offset = offset, limit = limit, scope = wide_scope)$data
            n <- nchar(right$Protein.Ids)
            stopifnot(all(n >= 2048L), all(substring(right$Protein.Ids, n - 2047L) == suffix))
            right$Protein.Ids <- substring(right$Protein.Ids, 1L, n - 2048L)
            stopifnot(identical(left, right))
            page_hashes[[i]] <- digest::digest(left, algo = "sha256")
        }
        releaseWorkflowTableReadScope(wide_scope)
        releaseWorkflowTableReadScope(scope)
        checkpoint$components$peptide_data <- wide
        withDIAStateScope(manager, function() commitDIAReferenceState(manager,
            "protein_peptide_filtered", checkpoint, list(), "PERF-106 wide screen"))
        wide_input <- primary
        wide_input$store <- manager$.context$store$directory
        wide_input$reference <- wide
        wide_input$head_sha256 <- digest::digest(file = file.path(wide_input$store,
            "dia-state-head.rds"), algo = "sha256")
        wide_input$payload_sha256 <- NULL
        wide_input$annotation_suffix_bytes <- 2048L
        closeDIAStateAdapter(manager)
        wide_input <- validate(wide_input, 622500L)
        stopifnot(identical(digest::digest(file = file.path(primary$store,
            primary$reference$payload$relative_path), algo = "sha256"), primary$payload_sha256))
        result <- list(schema = 1L, screen = list(primary = primary, wide = wide_input),
            full = full, build = workflowExecutionBuildIdentity(), library = normalizePath(library_path),
            preparation_page_hashes = page_hashes, source_screen_unchanged = TRUE,
            protocol_sha256 = digest::digest(file = file.path(repo,
                "tickets/openai/verification/perf106-repair-protocol-v2.json"), algo = "sha256"))
        saveRDS(result, output)
        receipt$status <- "complete"
        receipt$output <- output
        receipt$output_sha256 <- digest::digest(file = output, algo = "sha256")
    }, error = function(error) {
        receipt$status <<- "failed"
        receipt$error <<- conditionMessage(error)
    })
    receipt$finished <- Sys.time()
    flush()
    if (receipt$status != "complete") stop(receipt$error)
    receipt
}

#' Attribute serialization allocations and preflight separate result verification.
#' @param library_path,run_dir,repo Installed build, guarded output and frozen source.
#' @param verification_only Authorized row-name verifier repeat; reuse allocation evidence.
#' @return Preserved diagnostic receipt; not a timed performance member.
#' @noRd
checkPerf106RepairAttribution <- function(library_path, run_dir, repo, verification_only = FALSE) {
    binding <- jsonlite::read_json(file.path(repo,
        "tickets/openai/verification/pilot-source-bindings.json"))
    root <- file.path(binding$output_root, "repair-v2")
    ledger <- file.path(root, if (verification_only) "row-name-verifier-attempt.rds" else
        "attribution-attempt.rds")
    stopifnot(!file.exists(ledger))
    receipt <- list(status = "starting", run_dir = run_dir, started = Sys.time(),
        performance_claim = FALSE)
    flush <- function() {
        saveRDS(receipt, ledger)
        saveRDS(receipt, file.path(run_dir, "perf106-repair-attribution.rds"))
    }
    flush()
    tryCatch({
        receipt$verification_build <- workflowExecutionBuildIdentity()
        receipt$verification_library <- normalizePath(library_path)
        if (verification_only) {
            previous <- readRDS(file.path(root, "attribution-attempt.rds"))
            stopifnot(identical(previous$status, "failed"),
                identical(previous$error,
                    "identical(oracle, receipt$verification$selection[fields]) is not TRUE"))
            receipt$previous_attribution_sha256 <- digest::digest(
                file = file.path(root, "attribution-attempt.rds"), algo = "sha256")
        } else {
        installed <- readRDS(file.path(root, "correction-3/summary/application.rds"))
        stopifnot(identical(normalizePath(library_path), normalizePath(installed$library)),
            identical(workflowExecutionBuildIdentity(), installed$build),
            isTRUE(capabilities("profmem")))
        store <- newWorkflowEntityStore(run_dir, "perf106-attribution", "DIA")
        on.exit(closeWorkflowEntityStore(store), add = TRUE)
        value <- rep(strrep("x", 2048L), 10000L)
        measure <- function(name, fn) {
            path <- file.path(run_dir, paste0(name, "-allocations.txt"))
            utils::Rprofmem(path, threshold = 1048576)
            start <- proc.time()[["elapsed"]]
            hash <- tryCatch(fn(), finally = utils::Rprofmem(NULL))
            elapsed <- proc.time()[["elapsed"]] - start
            bytes <- suppressWarnings(as.numeric(sub(" .*", "", readLines(path))))
            list(hash = hash, seconds = elapsed,
                largest_R_allocation_bytes = max(c(0, bytes), na.rm = TRUE))
        }
        receipt$old_hash <- measure("whole-column", function()
            digest::digest(value, algo = "sha256", serializeVersion = 2L))
        receipt$streamed_hash <- measure("streamed-column", function()
            hashWorkflowSerializedValue(value, store))
        stopifnot(identical(receipt$old_hash$hash, receipt$streamed_hash$hash),
            receipt$old_hash$largest_R_allocation_bytes > 16 * 1024^2,
            receipt$streamed_hash$largest_R_allocation_bytes < 4 * 1024^2)
        receipt$allocation_scope <- paste("Parent R allocations only; child hashing processes",
            "are included in the enclosing service guard, not Rprofmem.",
            "This reproduces the serialization mechanism, not the historical OOM allocator.")
        }
        receipt$verification <- list()
        for (arm in c("native", "physical", "selection")) {
            member_dir <- file.path(root, "preflight-1", arm)
            path <- file.path(member_dir, "member-result.rds")
            before <- digest::digest(file = path, algo = "sha256")
            member <- readRDS(path)
            stopifnot(identical(member$status, "complete"), isTRUE(member$small_science_exact))
            project <- member$project_root
            owner_identity <- readSavedStoreIdentity(file.path(project, ".workflow-entities/store.json"))
            manager <- newDIAStateAdapter(newWorkflowEntityStore(project,
                owner_identity$project_id, owner_identity$workflow_id), reopen = TRUE)
            reference <- getDIAStateReference(manager)
            closeDIAStateAdapter(manager)
            receipt$verification[[arm]] <- verifyPerf106Result(project, reference)
            if (arm == "native") {
                wrong <- reference
                wrong$entity_id <- paste0("wrong-", reference$entity_id)
                error <- tryCatch({ verifyPerf106Result(project, wrong); NULL }, error = identity)
                stopifnot(inherits(error, "error"))
                receipt$verifier_rejects_invalid_reference <- list(
                    class = class(error), message = conditionMessage(error))
            }
            stopifnot(identical(digest::digest(file = path, algo = "sha256"), before))
        }
        fields <- c("rows", "samples", "columns", "page_rows", "page_hashes",
            "design_sha256", "summary_sha256", "support_sha256", "required_outputs_complete")
        oracle <- receipt$verification$native[fields]
        stopifnot(identical(oracle, receipt$verification$physical[fields]),
            identical(oracle, receipt$verification$selection[fields]))
        receipt$status <- "complete"
    }, error = function(error) {
        receipt$status <<- "failed"
        receipt$error <<- conditionMessage(error)
    })
    receipt$finished <- Sys.time()
    flush()
    if (receipt$status != "complete") stop(receipt$error)
    receipt
}

#' Execute a single budgeted installed qualification member.
#' @param library_path,run_dir,repo Installed package, owned member and source snapshot.
#' @param selection Fixed protocol member identifier.
#' @return Durable member receipt; never automatic promotion.
#' @noRd
runPerf106Qualification <- function(library_path, run_dir, repo, selection) {
    if (identical(selection, "perf106:V2-J"))
        return(checkPerf106ProfileInstrumentation(run_dir, repo))
    if (selection %in% c("perf106:V2-D1", "perf106:V2-D2"))
        return(verifyPerf106RepairP14(run_dir, repo,
            if (selection == "perf106:V2-D1") 21L else 22L, stage = "p14-diagnostic"))
    if (grepl("^perf106:V2-Q2[1-4]$", selection))
        return(verifyPerf106RepairP14(run_dir, repo,
            as.integer(sub("perf106:V2-Q", "", selection, fixed = TRUE))))
    if (grepl("^perf106:V2-Q([1-9]|1[0-9]|20)$", selection))
        return(verifyPerf106RepairScreen(library_path, run_dir, repo,
            as.integer(sub("perf106:V2-Q", "", selection, fixed = TRUE))))
    if (identical(selection, "perf106:V2-I"))
        return(preparePerf106RepairInputs(library_path, run_dir, repo))
    if (identical(selection, "perf106:V2-V"))
        return(checkPerf106RepairAttribution(library_path, run_dir, repo, verification_only = TRUE))
    if (identical(selection, "perf106:V2-A"))
        return(checkPerf106RepairAttribution(library_path, run_dir, repo))
    member <- sub("perf106:", "", selection, fixed = TRUE)
    binding <- jsonlite::read_json(file.path(repo,
        "tickets/openai/verification/pilot-source-bindings.json"))
    build <- jsonlite::read_json(file.path(repo,
        "docs/qa/workflow-execution/perf-106/candidate-build.json"))
    root <- binding$output_root
    stopifnot(identical(normalizePath(dirname(run_dir)), normalizePath(root)))
    ledger <- file.path(root, paste0("PERF106-", member, "-attempt.rds"))
    if (member == "S1" && file.exists(ledger)) {
        previous <- readRDS(ledger)
        # One replacement explicitly authorized by the user's 2026-10-04
        # `continue` after the completed operation's verifier-scope failure.
        stopifnot(identical(previous$status, "failed"),
            grepl('dbQuoteIdentifier', previous$error, fixed = TRUE),
            grepl('conn = "NULL"', previous$error, fixed = TRUE),
            identical(basename(previous$run_dir), "browser-smoke.XmL13i"))
        ledger <- file.path(root, "PERF106-S1-replacement-1-attempt.rds")
    }
    if (member == "R1" && file.exists(ledger)) {
        previous <- readRDS(ledger)
        # One replacement explicitly authorized on 2026-10-04 after the
        # accepted predecessor's pre-rollup capability-reporter failure.
        stopifnot(identical(previous$status, "failed"),
            identical(previous$p14_role, "control"),
            grepl('could not find function "workflowFastWriterInfo"',
                previous$error, fixed = TRUE),
            identical(basename(previous$run_dir), "browser-smoke.1qTSN4"))
        ledger <- file.path(root, "PERF106-R1-replacement-1-attempt.rds")
    }
    if (member == "R1" && file.exists(ledger)) {
        previous <- readRDS(ledger)
        # One final replacement explicitly authorized on 2026-10-04 after the
        # pre-rollup pkgload-only helper failure.
        stopifnot(identical(previous$status, "failed"),
            identical(previous$p14_role, "control"),
            grepl("No packages loaded with pkgload", previous$error, fixed = TRUE),
            identical(basename(previous$run_dir), "browser-smoke.5YuJR5"))
        ledger <- file.path(root, "PERF106-R1-final-replacement-attempt.rds")
    }
    if (member == "L0" && file.exists(ledger)) {
        previous <- readRDS(ledger)
        # One replacement explicitly authorized on 2026-10-04 after the
        # pre-data service-limit dispatch failure. Never replace its receipt.
        stopifnot(identical(previous$status, "failed_before_preparation"),
            identical(previous$preparation_started, FALSE),
            identical(previous$source_data_read, FALSE),
            identical(basename(previous$run_dir), "browser-smoke.CZ2wrF"))
        ledger <- file.path(root, "PERF106-L0-replacement-1-attempt.rds")
    }
    wide_recovery <- FALSE
    if (member == "L0" && file.exists(ledger)) {
        previous <- readRDS(ledger)
        # One wide-only recovery explicitly authorized on 2026-10-04. It is
        # eligible only after the retained 12 GiB generator failure and may not
        # regenerate or mutate the completed primary input.
        stopifnot(identical(previous$status, "failed"),
            grepl("11.9 GiB/12.0 GiB used", previous$error, fixed = TRUE),
            identical(basename(previous$run_dir), "browser-smoke.Xu0y9g"))
        ledger <- file.path(root, "PERF106-L0-wide-recovery-1-attempt.rds")
        wide_recovery <- TRUE
    }
    if (wide_recovery && file.exists(ledger)) {
        previous <- readRDS(ledger)
        # One final wide-only recovery explicitly authorized by the user's
        # 2026-10-04 `continue`, after the pre-input test-fixture failure.
        stopifnot(identical(previous$status, "failed"),
            identical(previous$recovery, "authorized-wide-only-after-Xu0y9g"),
            grepl('inherits(error, "error") is not TRUE', previous$error, fixed = TRUE),
            identical(basename(previous$run_dir), "browser-smoke.PQ7ihX"))
        ledger <- file.path(root, "PERF106-L0-wide-final-recovery-attempt.rds")
        receipt_recovery <- "authorized-final-wide-after-PQ7ihX"
    } else receipt_recovery <- if (wide_recovery)
        "authorized-wide-only-after-Xu0y9g" else NULL
    if (file.exists(ledger)) stop("PERF-106 member allowance already consumed")
    receipt <- list(status = "starting", member = member, run_dir = run_dir,
        started = Sys.time(), source = readRDS(file.path(run_dir, "source-manifest.rds")),
        promotion_authorized = FALSE, phases = list(), spill_attempts = NA_integer_)
    path <- file.path(run_dir, "perf106-result.rds")
    flush <- function() {
        for (destination in c(path, ledger)) {
            temporary <- paste0(destination, ".writing")
            saveRDS(receipt, temporary)
            if (!file.rename(temporary, destination)) stop("cannot publish member journal")
        }
    }
    flush()
    phase <- function(label, fn) {
        start <- proc.time()[["elapsed"]]
        on.exit({
            receipt$phases[[label]] <<- list(seconds = proc.time()[["elapsed"]] - start)
            flush()
        })
        fn()
    }
    tryCatch({
        phase("bindings", function() {
            runtime_binding <- list(path = build$library, files = build$package_files)
            if (grepl("^R[1-4]$", member)) {
                p14 <- readRDS(file.path(repo,
                    "docs/qa/workflow-execution/perf-106/p14-regression-bindings.rds"))
                role <- p14$order[[member]]
                runtime_binding <- p14[[role]]
                receipt$p14_role <<- role
                receipt$p14_binding_sha256 <<- digest::digest(file = file.path(repo,
                    "docs/qa/workflow-execution/perf-106/p14-regression-bindings.rds"),
                    algo = "sha256")
            }
            stopifnot(identical(normalizePath(library_path), normalizePath(runtime_binding$path)))
            for (name in names(runtime_binding$files)) stopifnot(identical(
                digest::digest(file = file.path(library_path, name), algo = "sha256"),
                runtime_binding$files[[name]]))
            for (item in list(binding$source_receipt, binding$design_authority,
                binding$native_source)) {
                source <- if (startsWith(item$path, "/")) item$path else file.path(repo, item$path)
                stopifnot(identical(digest::digest(file = source, algo = "sha256"), item$sha256))
            }
            writer <- if (exists("workflowFastWriterInfo", mode = "function", inherits = TRUE))
                workflowFastWriterInfo() else list(status = "accepted-hash-bound-predecessor",
                    fastduckdb = as.character(utils::packageVersion("fastduckdb")))
            receipt$runtime <- list(session = utils::sessionInfo(),
                exact_sum = diaRollupExtensionInfo(), writer = writer)
        })
        sys.source(file.path(repo, "tests/testthat/helper-workflow-bulk-rollup.R"),
            envir = environment())
        if (grepl("^S[1-6]$", member))
            sys.source(file.path(repo, "tools/profiling/workflow_perf101.R"),
                envir = environment())
        if (member == "L0") {
            receipt$inputs <- if (wide_recovery) {
                receipt$recovery <- receipt_recovery
                phase("wide_writer_lifecycle", checkPerf106WideWriterLifecycle)
                phase("prepare_and_check_wide", function()
                    recoverPerf106Wide(binding, run_dir, repo, copyDIAReplayCheckpoint))
            } else phase("prepare_and_check_inputs", function()
                preparePerf106Inputs(binding, run_dir, repo, copyDIAReplayCheckpoint))
        } else if (member == "L1") {
            receipt$lifecycle <- phase("installed_lifecycle", function()
                checkPerf106InstalledLifecycle(run_dir, repo))
        } else if (grepl("^S[1-6]$", member)) {
            receipt$screen <- runPerf106ScreenMember(member, binding, run_dir,
                copyDIAReplayCheckpoint, computePerf101Native)
        } else if (grepl("^F[1-8]$", member)) {
            receipt$confirmation <- runPerf106ConfirmationMember(member, run_dir, repo,
                copyDIAReplayCheckpoint)
        } else if (grepl("^R[1-4]$", member)) {
            receipt$p14_regression <- runPerf106P14Member(member, run_dir, repo,
                copyDIAReplayCheckpoint, fingerprintDIAWriterOutput)
        } else {
            stop("qualification member implementation is not yet frozen")
        }
        receipt$status <- "complete"
        receipt$finished <- Sys.time()
        flush()
        receipt
    }, error = function(error) {
        receipt$status <<- "failed"
        receipt$error <<- conditionMessage(error)
        receipt$finished <<- Sys.time()
        flush()
        stop(error)
    })
}

#' Run one accepted-P14 versus integrated-candidate regression member.
#' @param member Fixed R1--R4 member.
#' @param run_dir,repo Fresh member directory and frozen source snapshot.
#' @param copy_checkpoint,fingerprint Existing qualified helpers.
#' @return Complete operation, exact science and non-regression evidence.
#' @noRd
runPerf106P14Member <- function(member, run_dir, repo, copy_checkpoint, fingerprint) {
    binding <- readRDS(file.path(repo,
        "docs/qa/workflow-execution/perf-106/p14-regression-bindings.rds"))
    role <- binding$order[[member]]
    stopifnot(role %in% c("control", "candidate"),
        identical(digest::digest(file = binding$accepted_receipt, algo = "sha256"),
            binding$accepted_receipt_sha256),
        identical(digest::digest(file = file.path(binding$source_store,
            "dia-state-head.rds"), algo = "sha256"), binding$source_head_sha256))
    owned <- file.path(run_dir, "operation")
    dir.create(owned)
    started_at <- Sys.time()
    started <- proc.time()[["elapsed"]]
    manager <- copy_checkpoint(binding$source_store, owned)
    manager$.context$reference_states <- TRUE
    selectDIAStateReference(manager, "protein_peptide_filtered")
    checkpoint <- readDIAStateMetadata(manager)$data
    stopifnot(identical(checkpoint$components$peptide_data$entity_id,
        binding$input_reference$entity_id))
    phase_times <- list()
    phase <- function(label, fn) {
        begin <- proc.time()[["elapsed"]]
        on.exit(phase_times[[label]] <<- proc.time()[["elapsed"]] - begin)
        fn()
    }
    project <- file.path(run_dir, "dia_baseline_smoke")
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
    compute_p14 <- function(store, input, persist = NULL) {
        stopifnot(!is.null(persist))
        computeDIARelationalRollup(store, input, persist, "D", NULL, phase,
            diaRollupExtensionInfo()$path, profile = "P14")
    }
    result <- testthat::with_mocked_bindings(
        withDIAStateScope(manager, function() {
            result <- phase("apply_including_read_kernel_write_commit",
                function() runDIAReferencePrecursorRollup(workflow))
            arguments <- list(data = result$filteringSummary,
                step_name = "4_precursor_rollup", omic_type = "proteomics",
                experiment_label = "dia_baseline_smoke", return_grid = TRUE,
                overwrite = TRUE, project_dirs = directories)
            phase("publication", function() publishDIAQvalueOutput(manager, result,
                arguments, updateProteinFiltering))
            result
        }),
        computeSavedWorkflowPrecursorRollup = compute_p14,
        .package = "MultiScholaR"
    )
    record <- readDIAStateMetadata(manager)
    files <- list.files(directories$proteomics$time_dir, recursive = TRUE, full.names = TRUE)
    figures <- files[grepl("4_precursor_rollup_.*[.](png|pdf)$", files) &
        !grepl("/figures/", files)]
    stopifnot(length(figures) == 12L, all(file.size(figures) > 0),
        manager$.context$store$open_resources == 0L)
    value <- list(reference = record$data$components$peptide_data,
        summary = record$data$slots$args$rollUpPrecursorToPeptide$rollup_summary,
        filtering = unclass(result$filteringSummary),
        figures = stats::setNames(vapply(figures, file.size, numeric(1)), basename(figures)),
        execution = manager$.context$last_rollup_execution)
    closeDIAStateAdapter(manager)
    operation_seconds <- proc.time()[["elapsed"]] - started
    boundary <- c(list(started_at = started_at, finished_at = Sys.time(),
        elapsed_seconds = operation_seconds), readPerf106CgroupBoundary())
    saveRDS(boundary, file.path(run_dir, "operation-boundary.rds"))
    hashes <- fingerprint(owned, value$reference)
    summary_json <- jsonlite::toJSON(value$summary, auto_unbox = TRUE,
        null = "null", digits = NA)
    accepted_json <- jsonlite::toJSON(binding$accepted_summary, auto_unbox = TRUE,
        null = "null", digits = NA)
    stopifnot(identical(unname(hashes), unname(binding$accepted_output_hashes)),
        identical(summary_json, accepted_json), length(value$figures) == 12L,
        all(value$figures > 0),
        identical(value$execution$exact_sum_sha256,
            binding$accepted_exact_sum_sha256),
        identical(digest::digest(file = file.path(binding$source_store,
            "dia-state-head.rds"), algo = "sha256"), binding$source_head_sha256))
    list(member = member, role = role, profile = "P14",
        operation_boundary = boundary, natural_phase_seconds = phase_times,
        input_rows = 24000000L, output_rows = 12000000L,
        output_pages = length(hashes), science_exact = TRUE,
        summary_exact = TRUE, required_figures = length(value$figures),
        staged_bytes = value$reference$payload$bytes,
        source_unchanged = TRUE, spill_attempts = NA_integer_,
        performance_claim = FALSE)
}

#' Read exact cgroup counters at the complete-operation boundary.
#' @return Peak/current/swap values, before independent result verification.
#' @noRd
readPerf106CgroupBoundary <- function() {
    line <- grep("^0::/", readLines("/proc/self/cgroup"), value = TRUE)
    stopifnot(length(line) == 1L)
    group <- substring(line, 4L)
    root <- file.path("/sys/fs/cgroup", sub("^/", "", group))
    read_number <- function(name) {
        value <- readLines(file.path(root, name), n = 1L)
        stopifnot(length(value) == 1L, grepl("^[0-9]+$", value))
        as.numeric(value)
    }
    list(memory_current_bytes = read_number("memory.current"),
        memory_peak_bytes = read_number("memory.peak"),
        swap_peak_bytes = read_number("memory.swap.peak"))
}

#' Execute one frozen full-size native/finalist confirmation member.
#' @param member Fixed F1--F8 member.
#' @param run_dir,repo Fresh member directory and frozen source snapshot.
#' @param copy_checkpoint Existing replay copy helper.
#' @return Complete-operation receipt and separately bounded verification.
#' @noRd
runPerf106ConfirmationMember <- function(member, run_dir, repo, copy_checkpoint) {
    plan <- list(F1 = c("primary", "native"), F2 = c("primary", "physical"),
        F3 = c("primary", "physical"), F4 = c("primary", "native"),
        F5 = c("wide", "native"), F6 = c("wide", "physical"),
        F7 = c("wide", "physical"), F8 = c("wide", "native"))
    population <- unname(plan[[member]][[1L]])
    mode <- unname(plan[[member]][[2L]])
    decision <- jsonlite::read_json(file.path(repo,
        "docs/qa/workflow-execution/perf-106/screen-decision.json"))
    stopifnot(identical(decision$status, "complete"),
        identical(decision$finalist, "physical"),
        mode == "native" || identical(mode, decision$finalist))
    bindings <- readRDS(file.path(repo,
        "docs/qa/workflow-execution/perf-106/confirmation-input-bindings.rds"))
    prepared <- bindings[[population]]
    stopifnot(prepared$reference$metadata$shape[[1L]] == 6225000L,
        length(prepared$samples) == 300L,
        identical(digest::digest(file = file.path(prepared$store, "dia-state-head.rds"),
            algo = "sha256"), prepared$head_sha256),
        identical(digest::digest(file = file.path(prepared$store,
            prepared$reference$payload$relative_path), algo = "sha256"),
            prepared$payload_sha256))
    owned <- file.path(run_dir, "operation")
    dir.create(owned)
    started_at <- Sys.time()
    started <- proc.time()[["elapsed"]]
    manager <- copy_checkpoint(prepared$store, owned)
    manager$.context$reference_states <- TRUE
    input_head <- getDIAStateReference(manager)
    value <- if (mode == "native") computePerf106Native(manager, input_head,
        prepared$parameters) else computeDIAReferenceSampleFilter(manager, input_head,
            prepared$parameters, output_mode = mode)
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- manager
    workflow$qc_params <- list(peptide_qc = list())
    product <- withDIAStateScope(manager, function() finishDIAReferencePeptideStep(workflow,
        value$checkpoint, value$state_name, value$description, value$config, value$section,
        value$counts,
        observed_inclusion_absence = value$summary$inclusion_entries_absent_from_data))
    publication <- file.path(run_dir, "publication")
    dir.create(publication)
    publishWorkflowPeptideResult(list(manager = manager), list(request_dir = publication),
        list(state_reference = product$stateReference,
            operation = "minimum_peptides_per_sample"))
    stopifnot(is.null(assertDIAExportReady(manager)),
        manager$.context$store$open_resources == 0L)
    result_head <- getDIAStateReference(manager)
    closeDIAStateAdapter(manager)
    operation_seconds <- proc.time()[["elapsed"]] - started
    boundary <- c(list(started_at = started_at, finished_at = Sys.time(),
        elapsed_seconds = operation_seconds), readPerf106CgroupBoundary())
    saveRDS(boundary, file.path(run_dir, "operation-boundary.rds"))
    figures <- list.files(publication, pattern = "[.](png|pdf)$",
        recursive = TRUE, full.names = TRUE)
    stopifnot(length(figures) == 12L, all(file.size(figures) > 0),
        identical(digest::digest(file = file.path(prepared$store, "dia-state-head.rds"),
            algo = "sha256"), prepared$head_sha256))
    verification <- verifyPerf106Result(owned, result_head)
    stopifnot(verification$rows == 6150000L, verification$samples == 225L)
    list(member = member, population = population, mode = mode,
        operation_boundary = boundary,
        boundary = paste("adoption, parameter resolution, validation, compute, output/index/summary,",
            "checkpoint, figures, publication and teardown"),
        application_cache = "fresh request-owned copy; no filter result/count/selection cache",
        os_cache = "uncontrolled; no cold-disk claim or arm-specific prewarm",
        output_reference = verification$reference,
        staged_bytes = verification$reference$payload$bytes,
        figures = length(figures), verification = verification,
        spill_attempts = NA_integer_, performance_claim = FALSE)
}

#' Run the unchanged native S4 filter and existing persistence boundary at full size.
#' @param manager,reference,parameters Existing operation inputs.
#' @return Existing compute-only checkpoint envelope.
#' @noRd
computePerf106Native <- function(manager, reference, parameters) {
    checkpoint <- readDIAStateMetadata(manager, reference)$data
    object <- readDIAExecutionComponents(manager, reference,
        fields = c("peptide_data", "design_matrix"))
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
    checkpoint$components$peptide_data <- saveDIATableComponent(
        manager$.context, object@peptide_data)
    checkpoint$dimensions <- list(rows = nrow(object@peptide_data),
        columns = names(object@peptide_data))
    checkpoint$dependencies <- collectDIANativeReferences(list(reference,
        checkpoint$dependencies, checkpoint$components, checkpoint$support,
        checkpoint$products))
    list(checkpoint = checkpoint, reference = checkpoint$components$peptide_data,
        summary = summary, dependencies = checkpoint$dependencies, counts = NULL,
        execution = list(strategy = "unchanged native S4 filter; single process"),
        state_name = "sample_filtered", description = "Applied sample quality filter",
        config = parameters, section = "filterMinNumPeptidesPerSample")
}

#' Execute one frozen complete-operation screen arm without profiling wrappers.
#' @param member Fixed S1--S6 member.
#' @param binding Frozen pilot source authority.
#' @param run_dir Fresh guarded member directory.
#' @param copy_checkpoint Existing replay copy helper.
#' @param native_compute Unchanged native control adapter.
#' @return Complete-operation and independent bounded verification receipt.
#' @noRd
runPerf106ScreenMember <- function(member, binding, run_dir, copy_checkpoint,
    native_compute) {
    modes <- c(S1 = "native", S2 = "physical", S3 = "selection",
        S4 = "selection", S5 = "physical", S6 = "native")
    mode <- unname(modes[[member]])
    prepared_path <- file.path(binding$output_root, "screen-input.rds")
    prepared <- readRDS(prepared_path)
    stopifnot(prepared$reference$metadata$shape[[1L]] == 622500L,
        length(prepared$samples) == 300L,
        identical(digest::digest(file = file.path(prepared$store, "dia-state-head.rds"),
            algo = "sha256"), prepared$head_sha256))
    owned <- file.path(run_dir, "operation")
    dir.create(owned)
    started <- proc.time()[["elapsed"]]
    manager <- copy_checkpoint(prepared$store, owned)
    manager$.context$reference_states <- TRUE
    input_head <- getDIAStateReference(manager)
    value <- if (mode == "native") native_compute(manager, input_head,
        prepared$parameters) else computeDIAReferenceSampleFilter(manager, input_head,
            prepared$parameters, output_mode = mode)
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- manager
    workflow$qc_params <- list(peptide_qc = list())
    product <- withDIAStateScope(manager, function() finishDIAReferencePeptideStep(workflow,
        value$checkpoint, value$state_name, value$description, value$config, value$section,
        value$counts,
        observed_inclusion_absence = value$summary$inclusion_entries_absent_from_data))
    publication <- file.path(run_dir, "publication")
    dir.create(publication)
    publishWorkflowPeptideResult(list(manager = manager), list(request_dir = publication),
        list(state_reference = product$stateReference,
            operation = "minimum_peptides_per_sample"))
    stopifnot(is.null(assertDIAExportReady(manager)),
        manager$.context$store$open_resources == 0L)
    result_head <- getDIAStateReference(manager)
    closeDIAStateAdapter(manager)
    operation_seconds <- proc.time()[["elapsed"]] - started
    figures <- list.files(publication, pattern = "[.](png|pdf)$",
        recursive = TRUE, full.names = TRUE)
    stopifnot(length(figures) == 12L, all(file.size(figures) > 0),
        identical(digest::digest(file = file.path(prepared$store, "dia-state-head.rds"),
            algo = "sha256"), prepared$head_sha256))
    verification <- verifyPerf106Result(owned, result_head)
    stopifnot(verification$rows == 615000L, verification$samples == 225L)
    list(mode = mode, operation_seconds = operation_seconds,
        boundary = paste("adoption, parameter resolution, validation, compute, output/index/summary,",
            "checkpoint, figures, publication and teardown"),
        application_cache = "fresh request-owned copy; no filter result/count/selection cache",
        os_cache = "uncontrolled; no cold-disk claim or arm-specific prewarm",
        output_reference = verification$reference, staged_bytes = verification$reference$payload$bytes,
        figures = length(figures), verification = verification,
        spill_attempts = NA_integer_, performance_claim = FALSE)
}

#' Verify one accepted sample-filter result through fixed logical pages.
#' @param owned Completed request-owned project root.
#' @param state_reference Accepted state reference.
#' @return Stable table, design, summary and support fingerprints.
#' @noRd
verifyPerf106Result <- function(owned, state_reference) {
    identity <- readSavedStoreIdentity(file.path(owned, ".workflow-entities", "store.json"))
    manager <- newDIAStateAdapter(newWorkflowEntityStore(owned,
        identity$project_id, identity$workflow_id), reopen = TRUE)
    on.exit(closeDIAStateAdapter(manager), add = TRUE)
    record <- readDIAStateMetadata(manager, state_reference)$data
    reference <- record$components$peptide_data
    scope <- newWorkflowTableReadScope(manager$.context$store, reference, reader = "fastduckdb")
    on.exit(releaseWorkflowTableReadScope(scope), add = TRUE, after = FALSE)
    hashes <- character()
    for (offset in seq.int(0L, max(0L, scope$metadata$rows - 1L), by = 10000L)) {
        page <- readSavedWorkflowTable(manager$.context$store, reference,
            scope$metadata$columns, offset = offset,
            limit = min(10000L, scope$metadata$rows - offset), scope = scope)$data
        hashes[[length(hashes) + 1L]] <- digest::digest(page, algo = "sha256")
    }
    design <- materializeDIATable(manager$.context, record$components$design_matrix)
    support <- lapply(record$support$filterMinNumPeptidesPerSample,
        function(reference) materializeDIATable(manager$.context, reference))
    list(reference = reference, rows = scope$metadata$rows,
        samples = nrow(design), columns = scope$metadata$columns,
        page_rows = 10000L, page_hashes = hashes,
        design_sha256 = digest::digest(design, algo = "sha256"),
        summary_sha256 = digest::digest(
            record$slots$args$filterMinNumPeptidesPerSample$filter_summary, algo = "sha256"),
        support_sha256 = digest::digest(support, algo = "sha256"),
        required_outputs_complete = is.null(assertDIAExportReady(manager)),
        resources_released_after_close = TRUE)
}

#' Verify installed candidate completion, publication recovery and saved summaries.
#' @param run_dir,repo Owned member directory and frozen source snapshot.
#' @return Per-candidate scientific and output evidence.
#' @noRd
checkPerf106InstalledLifecycle <- function(run_dir, repo) {
    sys.source(file.path(repo, "tests/testthat/helper-workflow-checkpoint-fixtures.R"),
        envir = environment())
    object <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small", repo))
    samples <- as.character(object@design_matrix[[object@sample_id]])
    parameters <- list(peptides_per_sample_cutoff = 2L, inclusion_list = samples[[1L]])
    native <- filterMinNumPeptidesPerSample(object,
        peptides_per_sample_cutoff = parameters$peptides_per_sample_cutoff,
        inclusion_list = parameters$inclusion_list, core_utilisation = NA)
    results <- list()
    for (mode in c("physical", "selection")) {
        owned <- file.path(run_dir, paste0("lifecycle-", mode))
        dir.create(owned)
        manager <- newDIAStateAdapter(newWorkflowEntityStore(owned,
            paste0("perf106-lifecycle-", mode), "DIA"), WorkflowState$new(audit_enabled = FALSE))
        on.exit(closeDIAStateAdapter(manager), add = TRUE)
        manager$.context$reference_states <- TRUE
        manager$saveState("protein_peptide_filtered", object, list(), "PERF-106 lifecycle")
        input_head <- getDIAStateReference(manager)
        input <- readDIAStateMetadata(manager, input_head)$data$components$peptide_data
        input_path <- file.path(manager$.context$store$directory, input$payload$relative_path)
        input_sha256 <- digest::digest(file = input_path, algo = "sha256")
        value <- computeDIAReferenceSampleFilter(manager, input_head, parameters,
            output_mode = mode)
        actual <- materializeDIATable(manager$.context, value$reference)
        stopifnot(identical(actual, native@peptide_data),
            identical(value$summary,
                native@args$filterMinNumPeptidesPerSample$filter_summary),
            identical(getDIAStateReference(manager), input_head))
        inflated <- value
        inflated$counts$bounded_handoff_padding <- raw(5 * 1024^2)
        count_reference <- saveDIALargeFilteringCountsReference(manager$.context$store,
            inflated, object@protein_id_column)
        stopifnot(isSavedWorkflowReference(count_reference))
        workflow <- new.env(parent = emptyenv())
        workflow$state_manager <- manager
        workflow$qc_params <- list(peptide_qc = list())
        product <- withDIAStateScope(manager, function() finishDIAReferencePeptideStep(workflow,
            value$checkpoint, value$state_name, value$description, value$config, value$section,
            count_reference,
            observed_inclusion_absence = value$summary$inclusion_entries_absent_from_data))
        record <- readDIAStateMetadata(manager, product$stateReference)$data
        stopifnot(identical(materializeDIATable(manager$.context,
            record$components$peptide_data), native@peptide_data),
            identical(materializeDIATable(manager$.context,
                record$components$design_matrix), native@design_matrix),
            any(vapply(record$dependencies, function(reference)
                isSavedWorkflowReference(reference) &&
                    identical(reference$entity_id, count_reference$entity_id), logical(1))))
        publication <- file.path(owned, "publication")
        dir.create(publication)
        owner <- list(manager = manager)
        failure <- NULL
        options(multischolar.workflow_execution_boundary = function(stage) {
            if (identical(stage, "publication")) stop("injected PERF-106 publication failure")
        })
        failure <- tryCatch({
            publishWorkflowPeptideResult(owner, list(request_dir = publication),
                list(state_reference = product$stateReference,
                    operation = "minimum_peptides_per_sample"))
            NULL
        }, error = identity)
        options(multischolar.workflow_execution_boundary = NULL)
        stopifnot(inherits(failure, "error"),
            grepl("injected PERF-106 publication failure", conditionMessage(failure), fixed = TRUE))
        pending <- tryCatch({ assertDIAExportReady(manager); FALSE },
            saved_entity_outputs_pending = function(error) TRUE)
        stopifnot(pending)
        accepted_before_retry <- getDIAStateReference(manager)
        publishWorkflowPeptideResult(owner, list(request_dir = publication),
            list(state_reference = product$stateReference,
                operation = "minimum_peptides_per_sample"))
        stopifnot(is.null(assertDIAExportReady(manager)),
            identical(getDIAStateReference(manager), accepted_before_retry),
            identical(digest::digest(file = input_path, algo = "sha256"), input_sha256))
        figures <- list.files(publication, pattern = "[.](png|pdf)$",
            recursive = TRUE, full.names = TRUE)
        stopifnot(length(figures) == 12L, all(file.size(figures) > 0),
            manager$.context$store$open_resources == 0L)
        results[[mode]] <- list(science_exact = TRUE, output_mode = mode,
            count_reference = count_reference, publication_failure_observed = TRUE,
            publication_retry_without_science = TRUE, figures = length(figures),
            input_unchanged = TRUE, resources_released = TRUE)
        closeDIAStateAdapter(manager)
    }
    results
}

#' Prove installed fast-writer completion and failure cleanup before wide output.
#' @return TRUE after owned files and unfinished temporaries behave as required.
#' @noRd
checkPerf106WideWriterLifecycle <- function() {
    root <- tempfile("perf106-writer-")
    dir.create(root)
    on.exit(unlink(root, recursive = TRUE), add = TRUE)
    prototype <- data.frame(value = character())
    metadata <- workflowSavedTableMetadata(prototype, ".rownames", 3L)
    completed <- file.path(root, "completed.parquet")
    writer <- openWorkflowTableWriter(completed, metadata, backend = "fast_r")
    appendWorkflowTableWriter(writer, data.frame(value = "ok"))
    finishWorkflowTableWriter(writer)
    stopifnot(file.exists(completed))
    failed <- file.path(root, "failed.parquet")
    writer <- openWorkflowTableWriter(failed, metadata, backend = "fast_r")
    temporary <- writer$write_path
    error <- tryCatch({
        appendWorkflowTableWriter(writer, data.frame(other = "wrong schema"))
        NULL
    }, error = identity)
    stopifnot(inherits(error, "error"))
    abortWorkflowTableWriter(writer)
    stopifnot(!file.exists(failed),
        is.null(temporary) || !file.exists(temporary),
        identical(list.files(root), "completed.parquet"))
    TRUE
}

#' Stream the accepted primary input into the fixed annotation-wide population.
#' @param binding Frozen retained input authority.
#' @param run_dir,repo Owned recovery directory and frozen source snapshot.
#' @param copy_checkpoint Existing checkpoint copy helper.
#' @return Primary/wide manifests after exact bounded-page comparison.
#' @noRd
recoverPerf106Wide <- function(binding, run_dir, repo, copy_checkpoint) {
    descriptor <- readRDS(file.path(repo,
        "docs/qa/workflow-execution/perf-106/L0-preserved-primary.rds"))
    primary_path <- file.path(descriptor$store, descriptor$table$payload$relative_path)
    primary_sha256 <- "f8b33c6c175be20e1bae80dc7c7ab0e14765420b8388600fc3978019bdf02ae1"
    stopifnot(file.exists(primary_path), descriptor$table$metadata$shape[[1L]] == 6225000L,
        identical(digest::digest(file = primary_path, algo = "sha256"), primary_sha256))
    owned <- file.path(run_dir, "wide")
    dir.create(owned)
    manager <- copy_checkpoint(descriptor$store, owned)
    on.exit(closeDIAStateAdapter(manager), add = TRUE)
    manager$.context$reference_states <- TRUE
    checkpoint <- readDIAStateMetadata(manager)$data
    primary <- checkpoint$components$peptide_data
    stopifnot(identical(primary$entity_id, descriptor$table$entity_id))
    object <- readDIAExecutionComponents(manager, fields = "design_matrix")
    sample <- object@sample_id
    samples <- as.character(object@design_matrix[[sample]])
    limits <- rep(c(40000L, 2000L, 1000L, 40000L), each = 75L)
    stopifnot(length(samples) == 300L, !anyNA(samples), !anyDuplicated(samples))
    scope <- newWorkflowTableReadScope(manager$.context$store, primary, reader = "fastduckdb")
    on.exit(releaseWorkflowTableReadScope(scope), add = TRUE, after = FALSE)
    metadata <- scope$metadata
    stopifnot(metadata$rows == 6225000L, "Protein.Ids" %in% metadata$columns,
        identical(metadata$column_metadata$Protein.Ids$typeof, "character"))
    identity <- newSavedWorkflowScientificIdentity(primary$entity_id,
        checkpoint$components$design_matrix$entity_id,
        digest::digest(list(samples, limits, "wide", 2048L), algo = "sha256"),
        "perf106-wide-bounded-v1", manager$.context$software)
    suffix <- strrep("x", 2048L)
    wide <- saveWorkflowTable(scope$store, NULL, "perf106-wide-recovery", identity,
        metadata$key_column, metadata = metadata, write_fn = function(path, restoration) {
            writer <- openWorkflowTableWriter(path, restoration, backend = "fast_r")
            on.exit(abortWorkflowTableWriter(writer), add = TRUE)
            appendWorkflowTableWriter(writer, readSavedWorkflowTablePrototype(scope))
            streamSavedWorkflowTable(scope, metadata$columns, function(data, offset, row_names) {
                data$Protein.Ids <- paste0(data$Protein.Ids, suffix)
                appendWorkflowTableWriter(writer, data)
            })
            finishWorkflowTableWriter(writer)
        })
    wide_scope <- newWorkflowTableReadScope(manager$.context$store, wide, reader = "fastduckdb")
    on.exit(releaseWorkflowTableReadScope(wide_scope), add = TRUE, after = FALSE)
    pages <- 0L
    fingerprints <- character()
    for (offset in seq.int(0L, metadata$rows - 1L, by = 10000L)) {
        limit <- min(10000L, metadata$rows - offset)
        left <- readSavedWorkflowTable(scope$store, primary, metadata$columns,
            offset = offset, limit = limit, reader = "fastduckdb", scope = scope)$data
        right <- readSavedWorkflowTable(scope$store, wide, metadata$columns,
            offset = offset, limit = limit, reader = "fastduckdb", scope = wide_scope)$data
        n <- nchar(right$Protein.Ids, type = "bytes")
        stopifnot(all(n >= 2048L),
            all(substr(right$Protein.Ids, n - 2047L, n) == suffix))
        right$Protein.Ids <- substr(right$Protein.Ids, 1L, n - 2048L)
        stopifnot(identical(left, right))
        pages <- pages + 1L
        fingerprints[[pages]] <- digest::digest(left, algo = "sha256")
    }
    stopifnot(pages == ceiling(metadata$rows / 10000L),
        identical(digest::digest(file = primary_path, algo = "sha256"), primary_sha256))
    releaseWorkflowTableReadScope(wide_scope)
    releaseWorkflowTableReadScope(scope)
    checkpoint$components$peptide_data <- wide
    checkpoint$dimensions <- list(rows = metadata$rows, columns = metadata$columns)
    withDIAStateScope(manager, function() commitDIAReferenceState(manager,
        "protein_peptide_filtered", checkpoint, list(), "PERF-106 wide"))
    wide_path <- file.path(manager$.context$store$directory, wide$payload$relative_path)
    result <- list(primary = list(store = descriptor$store, reference = descriptor$table,
        samples = samples, limits = limits, synthetic = TRUE, annotation_suffix_bytes = 0L,
        parameters = list(peptides_per_sample_cutoff = 3000L,
            inclusion_list = samples[76:150]), payload_sha256 = primary_sha256),
        wide = list(store = manager$.context$store$directory, reference = wide,
            samples = samples, limits = limits, synthetic = TRUE, annotation_suffix_bytes = 2048L,
            parameters = list(peptides_per_sample_cutoff = 3000L,
                inclusion_list = samples[76:150]), payload_sha256 =
                digest::digest(file = wide_path, algo = "sha256")),
        verification = list(rows = metadata$rows, pages = pages,
            page_fingerprints = fingerprints, values_types_order_exact_after_suffix_removal = TRUE,
            primary_unchanged = TRUE))
    saveRDS(result, file.path(binding$output_root, "PERF106-inputs.rds"))
    closeDIAStateAdapter(manager)
    result
}

#' Prepare the fixed full primary and annotation-wide physical populations once.
#' @param binding Frozen retained P14 input authority.
#' @param run_dir,repo Owned preparation directory and frozen source snapshot.
#' @param copy_checkpoint Existing copyDIAReplayCheckpoint helper.
#' @return Input manifests including checks and immutable payload bindings.
#' @noRd
preparePerf106Inputs <- function(binding, run_dir, repo, copy_checkpoint) {
    stopifnot(identical(digest::digest(file = binding$peptide_source$path, algo = "sha256"),
        binding$peptide_source$expected_sha256))
    source_store <- dirname(dirname(dirname(binding$peptide_source$path)))
    source_head <- file.path(source_store, "dia-state-head.rds")
    before <- digest::digest(file = source_head, algo = "sha256")
    outputs <- list()
    for (population in c("primary", "wide")) {
        owned <- file.path(run_dir, population)
        dir.create(owned)
        manager <- copy_checkpoint(source_store, owned)
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
        design <- predecessor$components$design_matrix
        stopifnot(length(samples) == 300L, !anyNA(samples), !anyDuplicated(samples),
            identical(design$entity_id, authority$cases$primary$reference$identity$design),
            identical(design$payload$sha256, checkpoint$components$design_matrix$payload$sha256))
        scope <- newWorkflowTableReadScope(manager$.context$store, checkpoint$components$peptide_data)
        on.exit(releaseWorkflowTableReadScope(scope), add = TRUE, after = FALSE)
        configureWorkflowBulkQueryScope(scope)
        con <- scope$connection
        q <- \(x) as.character(DBI::dbQuoteIdentifier(con, x))
        qs <- \(x) as.character(DBI::dbQuoteString(con, x))
        input <- paste0("read_parquet(", qs(scope$path), ")")
        counts <- collectWorkflowBulkMetadata(scope, paste0("SELECT ", q(sample),
            " AS sample, count(*) AS n, count(DISTINCT (", q(object@protein_id_column),
            ",", q(object@peptide_sequence_column), ")) AS identities FROM ", input,
            " GROUP BY ", q(sample)))
        stopifnot(nrow(counts) == 300L, all(counts$n == 40000), all(counts$n == counts$identities))
        limits <- rep(c(40000L, 2000L, 1000L, 40000L), each = 75L)
        rules <- paste0("WHEN ", qs(samples), " THEN ", limits, collapse = " ")
        fields <- workflowTableStorageColumns(scope$metadata, scope$metadata$columns)
        projections <- q(fields)
        if (population == "wide") {
            stopifnot("Protein.Ids" %in% fields,
                identical(scope$metadata$column_metadata$Protein.Ids$typeof, "character"))
            # Injective suffix expansion of each original annotation identity.
            # Keys/quantities are unchanged; this is synthetic annotation stress.
            projections[fields == "Protein.Ids"] <-
                '"Protein.Ids" || repeat(\'x\', 2048) AS "Protein.Ids"'
        }
        key <- if (scope$metadata$key_column == ".rownames" &&
            isTRUE(scope$metadata$automatic_row_names)) paste0(
            ', CAST(row_number() OVER (ORDER BY ".multischolar_order") AS VARCHAR)',
            ' AS ".multischolar_row_key"') else ""
        sql <- paste0("SELECT ", paste(projections, collapse = ","), key,
            ', CAST(row_number() OVER (ORDER BY ".multischolar_order") AS INTEGER)',
            ' AS ".multischolar_order" FROM (SELECT *, row_number() OVER (PARTITION BY ',
            q(sample), ' ORDER BY ".multischolar_order") AS perf_rank FROM ', input,
            ") WHERE perf_rank <= CASE ", q(sample), " ", rules, " ELSE 0 END",
            ' ORDER BY ".multischolar_order"')
        metadata <- scope$metadata
        metadata$rows <- 6225000L
        identity <- newSavedWorkflowScientificIdentity(binding$peptide_source$entity_id,
            checkpoint$components$design_matrix$entity_id,
            digest::digest(list(samples, limits, population), algo = "sha256"),
            paste0("perf106-", population, "-v1"), manager$.context$software)
        selected <- saveWorkflowTableQuery(scope, sql, metadata, paste0("perf106-", population), identity)
        releaseWorkflowTableReadScope(scope)
        check <- newWorkflowTableReadScope(manager$.context$store, selected)
        on.exit(releaseWorkflowTableReadScope(check), add = TRUE, after = FALSE)
        configureWorkflowBulkQueryScope(check)
        con <- check$connection
        selected_table <- paste0("read_parquet(", qs(check$path), ")")
        observed <- collectWorkflowBulkMetadata(check, paste0("SELECT ", q(sample),
            " AS sample, count(*) AS n FROM ", selected_table, " GROUP BY ", q(sample)))
        stopifnot(identical(as.integer(observed$n[match(samples, observed$sample)]), limits),
            sum(observed$n) == 6225000, sum(limits[limits >= 3000L | samples %in% samples[76:150]]) == 6150000)
        releaseWorkflowTableReadScope(check)
        checkpoint$components$peptide_data <- selected
        checkpoint$dimensions <- list(rows = metadata$rows, columns = metadata$columns)
        withDIAStateScope(manager, function() commitDIAReferenceState(manager,
            "protein_peptide_filtered", checkpoint, list(), paste("PERF-106", population)))
        outputs[[population]] <- list(store = manager$.context$store$directory, reference = selected,
            samples = samples, limits = limits, synthetic = TRUE, annotation_suffix_bytes =
                if (population == "wide") 2048L else 0L,
            parameters = list(peptides_per_sample_cutoff = 3000L, inclusion_list = samples[76:150]),
            head_sha256 = digest::digest(file = file.path(manager$.context$store$directory,
                "dia-state-head.rds"), algo = "sha256"))
        closeDIAStateAdapter(manager)
    }
    stopifnot(identical(digest::digest(file = source_head, algo = "sha256"), before))
    saveRDS(outputs, file.path(binding$output_root, "PERF106-inputs.rds"))
    outputs
}
