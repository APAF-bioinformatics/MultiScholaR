# PERF-108 runtime comparison bindings and independent verification.
# Execution stays in the existing installed action runner and coordinator.

#' Expand the approved fresh-control cache qualification schedule.
#' @param bindings Frozen builds, scientific populations and output root.
#' @param stage One of cache-screen, cache-rollup or cache-filters.
#' @return Fixed cases; a filter schedule requires a recorded eligible rollup choice.
#' @noRd
perf108CacheCases <- function(bindings, stage) {
    stopifnot(stage %in% c("cache-screen", "cache-rollup", "cache-filters", "closure-filters"))
    closure <- stage == "closure-filters"
    screen <- stage == "cache-screen"
    populations <- if (stage %in% c("cache-filters", "closure-filters")) c("primary", "wide") else "p14"
    selected <- if (closure) "current" else if (stage == "cache-filters") readRDS(file.path(bindings$root,
        "cache-rollup/decision.rds"))$selected else NULL
    if (stage == "cache-filters") stopifnot(length(selected) == 1L, selected %in% c("current", "candidate"))
    order <- if (screen) c("current", "candidate", "candidate", "current") else
        if (stage == "cache-rollup") c("accepted", "current", "candidate", "candidate", "current", "accepted") else
            c("accepted", selected, selected, "accepted")
    cases <- list()
    for (population in populations) {
        repetitions <- c(accepted = 0L, current = 0L, candidate = 0L)
        input <- bindings[[if (screen) "screen" else "full"]][[population]]
        for (role in order) {
            repetitions[[role]] <- repetitions[[role]] + 1L
            rep <- repetitions[[role]]
            cases[[length(cases) + 1L]] <- list(id = paste(population, role, rep, sep = "-"),
                variant = if (population == "p14") "p14" else "physical",
                role = role, control = role == "accepted", fault = "none", rep = rep,
                population = population, library = bindings[[role]]$library,
                expected_build = bindings[[role]]$build, profile_worker = FALSE,
                binding = digest::digest(list(input, bindings[[role]], stage, rep,
                    protocol = bindings$protocol_sha256), algo = "sha256"),
                input = input, verification_slot = length(cases) + 1L)
        }
    }
    cases
}

#' Decide from completed exact-output verifications and paired complete-operation measurements.
#' @param receipt,cases Current stage outcomes and its fixed cases.
#' @return Explicit scientific validity, comparison and selected role; never default activation.
#' @noRd
perf108CacheDecision <- function(receipt, cases) {
    eligible <- function(id) {
        item <- receipt$cases[[id]]
        op <- item$operation
        isTRUE(item$passed) && isTRUE(item$verifier_completion$accepted) &&
            identical(op$status, "complete") && isTRUE(op$measurement_valid) &&
            isTRUE(op$lifecycle_passed) && identical(op$swap_bytes, 0) &&
            is.finite(op$peak_bytes) && op$peak_bytes > 0 && op$peak_bytes < 15.5 * 1024^3 &&
            is.finite(op$elapsed_seconds) && op$elapsed_seconds > 0 &&
            op$elapsed_seconds <= if (receipt$stage == "cache-screen") 300 else 1800
    }
    compare <- function(population, role, reference) {
        pairs <- lapply(1:2, function(rep) {
            ids <- paste(population, c(role, reference), rep, sep = "-")
            valid <- all(vapply(ids, eligible, logical(1)))
            if (!valid) return(list(valid = FALSE, within = FALSE))
            x <- receipt$cases[[ids[[1L]]]]$operation
            y <- receipt$cases[[ids[[2L]]]]$operation
            peak <- x$peak_bytes / y$peak_bytes
            elapsed <- x$elapsed_seconds / y$elapsed_seconds
            list(valid = TRUE, within = peak <= 1.05 && elapsed <= 1.05,
                peak_ratio = peak, elapsed_ratio = elapsed,
                seconds = x$elapsed_seconds, reference_seconds = y$elapsed_seconds)
        })
        valid <- all(vapply(pairs, `[[`, logical(1), "valid"))
        list(valid = valid, within = all(vapply(pairs, `[[`, logical(1), "within")),
            median_elapsed_ratio = if (valid) stats::median(vapply(pairs, `[[`, numeric(1), "seconds")) /
                stats::median(vapply(pairs, `[[`, numeric(1), "reference_seconds")) else NULL, pairs = pairs)
    }
    valid <- all(vapply(vapply(cases, `[[`, character(1), "id"), eligible, logical(1)))
    result <- list(status = "held", science_outputs_and_measurements_valid = valid,
        selected = NULL, promotion_authorized = FALSE)
    if (receipt$stage == "cache-screen") {
        result$status <- if (valid) "screen_valid" else "screen_invalid"
        return(result)
    }
    if (receipt$stage == "cache-rollup") {
        result$current <- compare("p14", "current", "accepted")
        result$candidate <- compare("p14", "candidate", "accepted")
        result$cache <- compare("p14", "candidate", "current")
        result$meaningful_cache_benefit <- isTRUE(result$cache$within) &&
            result$cache$median_elapsed_ratio <= 0.90
        if (result$candidate$within && (!result$current$within || result$meaningful_cache_benefit))
            result$selected <- "candidate" else if (result$current$within) result$selected <- "current"
        if (!is.null(result$selected)) result$status <- "rollup_route_eligible"
    } else {
        selected <- unique(vapply(Filter(\(x) x$role != "accepted", cases), `[[`, character(1), "role"))
        stopifnot(length(selected) == 1L)
        result$comparisons <- stats::setNames(lapply(c("primary", "wide"),
            \(p) compare(p, selected, "accepted")), c("primary", "wide"))
        if (all(vapply(result$comparisons, `[[`, logical(1), "within"))) {
            result$status <- "runtime_non_regression_passed"
            result$selected <- selected
        }
    }
    result
}

#' Record phase boundaries and complete COPY query profiles in an untimed worker.
#' @param directory Existing request directory owned by the diagnostic member.
#' @return Finalizer restoring the prior option and closing the event log.
#' @noRd
startPerf108Attribution <- function(directory) {
    events <- file(file.path(directory, "phase-events.tsv"), "wt")
    writeLines("uptime\tid\tphase\tevent", events)
    clock <- function() as.numeric(strsplit(readLines("/proc/uptime", n = 1L), " ")[[1L]][[1L]])
    serial <- 0L
    emit <- function(id, label, event) {
        writeLines(paste(clock(), id, label, event, sep = "\t"), events)
        flush(events)
    }
    phase <- function(label, fn) {
        serial <<- serial + 1L
        id <- serial
        status <- "failed"
        emit(id, label, "start")
        on.exit(emit(id, label, status), add = TRUE)
        if (identical(label, "ordered_result_stage")) {
            frame <- environment(fn)
            con <- get("con", frame, inherits = TRUE)
            query <- get("query", frame, inherits = TRUE)
            plan <- workflowQueryCollect(con, paste0("EXPLAIN (FORMAT JSON) ", query))$explain_value
            writeLines(plan, file.path(directory, paste0("complete-plan-", id, ".json")))
            # Re-enabling profiling can flush the preceding disable statement.
            # Keep the engine's mutable destination separate from sealed COPY evidence.
            output <- file.path(directory, "active-query-profile.json")
            workflowQueryExecute(con, "SET profiling_coverage = 'ALL'")
            workflowQueryExecute(con, "SET enable_profiling = 'json'")
            on.exit(workflowQueryExecute(con, "PRAGMA disable_profiling"), add = TRUE)
            workflowQueryExecute(con, paste0("SET profiling_output = ", workflowQueryQuoteString(con, output)))
        }
        value <- fn()
        if (identical(label, "ordered_result_stage")) {
            stopifnot(startsWith(jsonlite::read_json(output)$query_name, "COPY"),
                file.copy(output, file.path(directory, paste0("query-profile-", id, ".json"))))
        }
        status <- "complete"
        value
    }
    previous <- options(multischolar.dia_rollup_phase = phase)
    emit(0L, "worker", "start")
    function() {
        emit(0L, "worker", "finished")
        options(previous)
        close(events)
    }
}

#' Add phase observation to the existing process-entry profiler.
#' @param source Frozen helper path; never a mutable scientific function override.
#' @param entry Fixed installed worker entry, or the bounded lifecycle fixture.
#' @return Worker expression for the existing diagnostic launcher.
#' @noRd
perf108AttributionExpression <- function(source,
    entry = "MultiScholaR:::workflowWorkerMain(commandArgs(TRUE))") {
    setup <- paste0('e <- new.env(parent=asNamespace("MultiScholaR")); sys.source(',
        encodeString(source, quote = '"'), ', e); ',
        'finish <- e$startPerf108Attribution(d); tryCatch({ ', entry,
        ' }, finally=finish())')
    perf106ProfileExpression(setup)
}

#' Expand the fixed old/new/new/old runtime comparison order.
#' @param bindings Frozen input/build/oracle bindings.
#' @param stage runtime-screen or runtime-full.
#' @return Fixed cases, never an adaptive candidate sweep.
#' @noRd
perf108Cases <- function(bindings, stage) {
    stopifnot(stage %in% c("runtime-screen", "runtime-full"))
    screen <- stage == "runtime-screen"
    populations <- if (screen) "p14" else c("p14", "primary", "wide")
    result <- list()
    for (population in populations) for (i in seq_len(4L)) {
        role <- c("control", "candidate", "candidate", "control")[[i]]
        rep <- c(1L, 1L, 2L, 2L)[[i]]
        input <- bindings[[if (screen) "screen" else "full"]][[population]]
        result[[length(result) + 1L]] <- list(id = paste(population, role, rep, sep = "-"),
            variant = if (population == "p14") "p14" else "physical",
            role = role, control = role == "control", fault = "none", rep = rep,
            population = population, library = bindings[[role]]$library,
            expected_build = bindings[[role]]$build,
            binding = digest::digest(list(input, bindings[[role]], stage), algo = "sha256"),
            input = input, verification_slot = length(result) + 1L + if (screen) 0L else 4L)
    }
    result
}

#' Expand only the corrected P14 candidate members; retained controls are reused.
#' @param bindings Frozen repaired-candidate and retained-control bindings.
#' @param stage runtime-repair-screen or runtime-repair-full.
#' @return Two fixed candidate cases, with no repeated control work.
#' @noRd
perf108RepairCases <- function(bindings, stage) {
    stopifnot(stage %in% c("runtime-repair-screen", "runtime-repair-screen-recheck",
        "runtime-repair-full", "runtime-staged-screen", "runtime-staged-full",
        "runtime-attribution-screen", "runtime-attribution-full",
        "runtime-batch-screen", "runtime-batch-full"))
    attribution <- startsWith(stage, "runtime-attribution-")
    batch <- startsWith(stage, "runtime-batch-")
    screen <- !stage %in% c("runtime-repair-full", "runtime-staged-full", "runtime-attribution-full", "runtime-batch-full")
    input <- bindings[[if (screen) "screen" else "full"]]$p14
    reps <- if (attribution || batch) 1L else if (stage == "runtime-repair-screen-recheck") 2L else seq_len(2L)
    lapply(reps, function(rep) list(
        id = paste("p14", "candidate", rep, sep = "-"),
        variant = "p14", role = "candidate", control = FALSE,
        fault = "none", rep = rep, population = "p14",
        profile_worker = attribution, attribute_worker = attribution,
        batch_count = if (batch) 6L else 3L,
        library = bindings$candidate$library,
        expected_build = bindings$candidate$build,
        binding = digest::digest(list(input, bindings$candidate, stage, rep), algo = "sha256"),
        input = input, verification_slot = rep
    ))
}

#' Verify one runtime member outside the measured application slice.
#' @param run_dir,repo Guarded verifier output and source snapshot.
#' @param slot Fixed screen/full slot, 1 through 16.
#' @return Source-bound verification receipt; no measured-operation work is repeated.
#' @noRd
verifyPerf108Runtime <- function(run_dir, repo, slot, canonical = FALSE, cache_stage = NULL,
    structural_stage = NULL) {
    stopifnot(slot %in% seq_len(16L))
    ordered <- !is.null(structural_stage) && startsWith(structural_stage, "ordered-")
    bindings <- readRDS(file.path(repo, "docs/qa/workflow-execution/perf-108",
        if (ordered) "ordered-groups/bindings.rds" else if (!is.null(structural_stage)) "structural/bindings.rds" else
            if (identical(cache_stage, "closure-filters")) "closeout/bindings.rds" else
            if (is.null(cache_stage)) "runtime-bindings.rds" else "cache-qualification/bindings.rds"))
    stage <- if (!is.null(structural_stage)) structural_stage else
        if (!is.null(cache_stage)) cache_stage else if (slot <= 4L) "runtime-screen" else "runtime-full"
    cases <- if (!is.null(structural_stage)) perf108StructuralCases(bindings, stage) else
        if (is.null(cache_stage)) perf108Cases(bindings, stage) else perf108CacheCases(bindings, stage)
    case <- cases[[if (!is.null(structural_stage) || !is.null(cache_stage) || slot <= 4L) slot else slot - 4L]]
    screen <- stage %in% c("runtime-screen", "cache-screen", "structural-screen", "ordered-screen")
    owned <- file.path(bindings$root, stage, case$id)
    verification_owned <- if (canonical) file.path(owned, "canonical-verification") else owned
    path <- file.path(verification_owned, "verification.rds")
    stopifnot(!file.exists(path))
    receipt <- list(status = "starting", run_dir = run_dir, member_id = case$id)
    on.exit({
        if (identical(receipt$status, "starting")) receipt$status <- "failed"
        saveRDS(receipt, path)
    }, add = TRUE)
    sys.source(file.path(repo, "tools/profiling/workflow_perf106.R"), environment())
    sys.source(file.path(repo, "tests/testthat/helper-workflow-bulk-rollup.R"), environment())
    operation_path <- file.path(owned, "operation-receipt.rds")
    receipt$operation_sha256 <- digest::digest(file = operation_path, algo = "sha256")
    operation <- readRDS(operation_path)
    member <- readRDS(file.path(owned, "member-result.rds"))
    stopifnot(identical(operation$status, "complete"), identical(member$status, "complete"),
        identical(operation$config$build, case$expected_build), identical(operation$binding, case$binding),
        identical(member$project_root, file.path(owned, "project")),
        identical(operation$member$parameters, case$input$parameters),
        identical(operation$member$input_entity_id, case$input$reference$entity_id),
        identical(digest::digest(file = file.path(case$input$store, "dia-state-head.rds"),
            algo = "sha256"), case$input$head_sha256))
    project <- member$project_root
    owner <- readSavedStoreIdentity(file.path(project, ".workflow-entities/store.json"))
    manager <- newDIAStateAdapter(newWorkflowEntityStore(project, owner$project_id, owner$workflow_id), reopen = TRUE)
    on.exit(closeDIAStateAdapter(manager), add = TRUE)
    checkpoint <- readDIAStateMetadata(manager)$data
    request <- file.path(project, ".workflow-entities/.workflow-export", operation$execution$request_id)
    worker <- readRDS(file.path(request, "worker-result.rds"))
    execution <- worker$result$execution
    if (case$population == "p14") {
        oracle <- bindings$p14_oracle
        stopifnot(identical(digest::digest(file = oracle$receipt, algo = "sha256"), oracle$receipt_sha256))
        reference <- checkpoint$components$peptide_data
        if (canonical) {
            receipt$original_physical_fingerprint_count <- length(fingerprintDIAWriterOutput(project, reference))
        }
        hashes <- fingerprintPerf108LogicalPages(project, reference)
        receipt$page_policy <- "fixed logical 10000-row pages; original golden hashes unchanged"
        if (canonical) {
            # Check the same logical page convention against the matched retained control.
            control_project <- file.path(bindings$root, stage,
                paste(case$population, "control", case$rep, sep = "-"), "project")
            control_owner <- readSavedStoreIdentity(file.path(control_project, ".workflow-entities/store.json"))
            control_manager <- newDIAStateAdapter(newWorkflowEntityStore(control_project,
                control_owner$project_id, control_owner$workflow_id), reopen = TRUE)
            on.exit(closeDIAStateAdapter(control_manager), add = TRUE)
            control_reference <- readDIAStateMetadata(control_manager)$data$components$peptide_data
            receipt$canonical_control_pages <- fingerprintPerf108LogicalPages(control_project, control_reference)
        }
        expected <- oracle[[if (screen) "screen_pages" else "full_pages"]]
        encode <- function(value) jsonlite::toJSON(value, auto_unbox = TRUE, null = "null", digits = NA)
        summary <- encode(checkpoint$slots$args$rollUpPrecursorToPeptide$rollup_summary)
        if (canonical) stopifnot(identical(unname(receipt$canonical_control_pages), unname(expected)))
        figures <- list.files(file.path(request, "required-output"),
            "4_precursor_rollup_.*[.](png|pdf)$", recursive = TRUE, full.names = TRUE)
        figures <- figures[!grepl("/figures/", figures)]
        stopifnot(identical(unname(hashes), unname(expected)),
            reference$metadata$shape[[1L]] == if (screen) 480000L else 12000000L,
            identical(execution$exact_sum_sha256, bindings[[case$role]]$exact_sum),
            identical(execution$profile$name, "P14"), identical(execution$engine_threads, 8L),
            identical(execution$writer_backend, "fast_r"), identical(execution$writer_threads, 8L),
            !is.null(execution$admission$worker_cap), length(figures) == 12L,
            all(file.size(figures) > 0), is.null(assertDIAExportReady(manager)))
        if (case$role %in% c("candidate", "current")) stopifnot(identical(execution$query_owner, "fastduckdb_c"))
        if (!is.null(cache_stage) || !is.null(structural_stage)) {
            stopifnot(identical(execution$profile$batches, 3L),
                identical(execution$batch_sizes, rep(if (screen) 4L else 100L, 3L)))
            if (case$role != "accepted") stopifnot(identical(execution$fastduckdb_threads, 8L))
        }
        if (ordered && case$role == "candidate") {
            stopifnot(identical(execution$aggregation_frontend, "ordered_native_groups"),
                identical(execution$ordered_numeric_threads, 1L),
                identical(execution$output_transport, "staged_or_R"), length(execution$ordered_reducer) == 3L,
                sum(vapply(execution$ordered_reducer, `[[`, numeric(1), "input_rows")) == if (screen) 960000 else 24000000,
                sum(vapply(execution$ordered_reducer, `[[`, numeric(1), "groups")) == if (screen) 480000 else 12000000,
                all(vapply(execution$ordered_reducer, `[[`, numeric(1), "max_active_groups") == 1),
                all(vapply(execution$ordered_reducer, `[[`, numeric(1), "max_group_state_bytes") <= 67108864))
        } else if (!is.null(structural_stage) && case$role == "candidate") stopifnot(
            identical(execution$output_transport, "owned_query_arrow_to_native_writer"),
            identical(execution$intermediate_bytes_written, 0))
        if (!screen) stopifnot(identical(summary, encode(oracle$full_summary)))
        receipt$summary <- summary
        receipt$pages <- hashes
        receipt$figures <- stats::setNames(file.size(figures), basename(figures))
    } else {
        oracle <- bindings$filter_oracles[[case$population]]
        stopifnot(identical(digest::digest(file = oracle$path, algo = "sha256"), oracle$sha256))
        verification <- verifyPerf106Result(project, getDIAStateReference(manager))
        fields <- c("rows", "samples", "columns", "page_rows", "page_hashes", "design_sha256",
            "summary_sha256", "support_sha256", "required_outputs_complete")
        stopifnot(identical(verification[fields], oracle$value[fields]),
            identical(execution$bulk_query$settings, list(threads = 8L, engine_limit = "12GiB", spill = FALSE)))
        receipt$verification <- verification
    }
    stopifnot(identical(digest::digest(file = operation_path, algo = "sha256"), receipt$operation_sha256))
    receipt$execution <- execution
    receipt$status <- "complete"
    receipt
}

#' Verify a corrected P14 candidate without repeating its retained control.
#' @param run_dir,repo Guarded verifier output and frozen source snapshot.
#' @param rep Candidate repetition, one or two.
#' @param full Use the full 300-sample input; otherwise use the retained screen.
#' @return Exact scientific/output verification bound to the corrected operation.
#' @noRd
verifyPerf108Repair <- function(run_dir, repo, rep, full = FALSE, recheck = FALSE,
    staged = FALSE, attribution = FALSE, batch = FALSE) {
    stopifnot(rep %in% 1:2, is.logical(full), length(full) == 1L)
    bindings <- readRDS(file.path(repo, "docs/qa/workflow-execution/perf-108",
        if (batch) "runtime-batch-bindings.rds" else if (attribution) "runtime-attribution-bindings.rds" else
            if (staged) "runtime-staged-bindings.rds" else "runtime-repair-bindings.rds"))
    stage <- if (batch) {
        if (full) "runtime-batch-full" else "runtime-batch-screen"
    } else if (attribution) {
        if (full) "runtime-attribution-full" else "runtime-attribution-screen"
    } else if (staged) {
        if (full) "runtime-staged-full" else "runtime-staged-screen"
    } else if (full) "runtime-repair-full" else if (recheck)
        "runtime-repair-screen-recheck" else "runtime-repair-screen"
    cases <- perf108RepairCases(bindings, stage)
    at <- match(rep, vapply(cases, `[[`, integer(1), "rep"))
    stopifnot(!is.na(at))
    case <- cases[[at]]
    owned <- file.path(bindings$root, stage, case$id)
    path <- file.path(owned, "verification.rds")
    stopifnot(!file.exists(path))
    receipt <- list(status = "starting", run_dir = run_dir, member_id = case$id)
    on.exit({
        if (identical(receipt$status, "starting")) receipt$status <- "failed"
        saveRDS(receipt, path)
    }, add = TRUE)
    sys.source(file.path(repo, "tools/profiling/workflow_perf106.R"), environment())
    operation_path <- file.path(owned, "operation-receipt.rds")
    receipt$operation_sha256 <- digest::digest(file = operation_path, algo = "sha256")
    operation <- readRDS(operation_path)
    member <- readRDS(file.path(owned, "member-result.rds"))
    stopifnot(identical(operation$status, "complete"), identical(member$status, "complete"),
        identical(operation$config$build, case$expected_build), identical(operation$binding, case$binding),
        identical(member$project_root, file.path(owned, "project")),
        identical(operation$member$input_entity_id, case$input$reference$entity_id),
        identical(digest::digest(file = file.path(case$input$store, "dia-state-head.rds"),
            algo = "sha256"), case$input$head_sha256))
    project <- member$project_root
    owner <- readSavedStoreIdentity(file.path(project, ".workflow-entities/store.json"))
    manager <- newDIAStateAdapter(newWorkflowEntityStore(project,
        owner$project_id, owner$workflow_id), reopen = TRUE)
    on.exit(closeDIAStateAdapter(manager), add = TRUE)
    checkpoint <- readDIAStateMetadata(manager)$data
    reference <- checkpoint$components$peptide_data
    hashes <- fingerprintPerf108LogicalPages(project, reference)
    expected <- bindings$p14_oracle[[if (full) "full_pages" else "screen_pages"]]
    request <- file.path(project, ".workflow-entities/.workflow-export", operation$execution$request_id)
    worker <- readRDS(file.path(request, "worker-result.rds"))
    execution <- worker$result$execution
    if (batch) stopifnot(identical(execution$profile$batches, 6L),
        identical(execution$batch_sizes, rep(if (full) 50L else 2L, 6L)))
    encode <- function(value) jsonlite::toJSON(value, auto_unbox = TRUE,
        null = "null", digits = NA)
    summary <- encode(checkpoint$slots$args$rollUpPrecursorToPeptide$rollup_summary)
    figures <- list.files(file.path(request, "required-output"),
        "4_precursor_rollup_.*[.](png|pdf)$", recursive = TRUE, full.names = TRUE)
    figures <- figures[!grepl("/figures/", figures)]
    stopifnot(identical(unname(hashes), unname(expected)),
        reference$metadata$shape[[1L]] == if (full) 12000000L else 480000L,
        identical(execution$exact_sum_sha256, bindings$candidate$exact_sum),
        identical(execution$query_owner, "fastduckdb_c"),
        identical(execution$profile$name, "P14"), identical(execution$engine_threads, 8L),
        identical(execution$writer_backend, "fast_r"), identical(execution$writer_threads, 8L),
        !is.null(execution$admission$worker_cap), length(figures) == 12L,
        all(file.size(figures) > 0), is.null(assertDIAExportReady(manager)))
    if (full) stopifnot(identical(summary, encode(bindings$p14_oracle$full_summary)))
    receipt$status <- "complete"
    receipt$summary <- summary
    receipt$pages <- hashes
    receipt$figures <- stats::setNames(file.size(figures), basename(figures))
    receipt$execution <- execution
    receipt
}

#' Compare corrected candidates with the retained control receipts.
#' @param receipt,cases,bindings Corrected stage, its cases and frozen bindings.
#' @return Scientific validity plus the unchanged full-size non-regression rule.
#' @noRd
perf108RepairDecision <- function(receipt, cases, bindings) {
    helpers <- new.env(parent = baseenv())
    sys.source(file.path(bindings$root, receipt$stage,
        "source/tools/profiling/workflow_perf106.R"), helpers)
    screen <- !receipt$stage %in% c("runtime-repair-full", "runtime-staged-full")
    control_stage <- if (screen) "runtime-screen" else "runtime-full"
    pairs <- list()
    valid <- TRUE
    regression <- TRUE
    for (rep in unique(vapply(cases, `[[`, integer(1), "rep"))) {
        control_id <- paste("p14", "control", rep, sep = "-")
        candidate_id <- paste("p14", "candidate", rep, sep = "-")
        control_dir <- file.path(bindings$control_root, control_stage, control_id)
        candidate_dir <- file.path(bindings$root, receipt$stage, candidate_id)
        control_stage_receipt <- readRDS(file.path(bindings$control_root,
            control_stage, "stage.rds"))
        control_case <- control_stage_receipt$cases[[control_id]]
        candidate_case <- receipt$cases[[candidate_id]]
        control_checked <- helpers$readPerf106Verification(control_dir,
            dirname(bindings$control_root), control_case$verifier_exit)
        candidate_checked <- helpers$readPerf106Verification(candidate_dir,
            dirname(bindings$root), candidate_case$verifier_exit)
        eligible <- function(operation, checked) isTRUE(checked$accepted) &&
            identical(operation$status, "complete") && isTRUE(operation$measurement_valid) &&
            isTRUE(operation$lifecycle_passed) && operation$peak_bytes < 15.5 * 1024^3 &&
            identical(operation$swap_bytes, 0) && operation$elapsed_seconds <= if (screen) 300 else 1800
        paired <- eligible(control_case$operation, control_checked) &&
            eligible(candidate_case$operation, candidate_checked) &&
            identical(control_checked$value$summary, candidate_checked$value$summary)
        ratios <- if (paired) list(
            peak = candidate_case$operation$peak_bytes / control_case$operation$peak_bytes,
            elapsed = candidate_case$operation$elapsed_seconds / control_case$operation$elapsed_seconds
        ) else NULL
        within <- paired && ratios$peak <= bindings$thresholds$peak &&
            ratios$elapsed <= bindings$thresholds$elapsed
        pairs[[paste0("p14-", rep)]] <- list(valid = paired,
            non_regression = within, ratios = ratios)
        valid <- valid && paired
        regression <- regression && within
    }
    list(status = if (!valid) "invalid_or_incomplete" else if (screen)
        "screen_valid" else if (regression) "non_regression_passed" else "non_regression_failed",
        science_outputs_and_measurements_valid = valid,
        non_regression = regression, pairs = pairs, promotion_authorized = FALSE,
        controls_reused = TRUE)
}

#' Apply the frozen runtime non-regression rule to verified member pairs.
#' @param receipt,cases,bindings Existing measurements, fixed cases and frozen policy.
#' @return Decision with separate scientific validity and measured non-regression.
#' @noRd
perf108Decision <- function(receipt, cases, bindings, verification_overrides = list()) {
    helpers <- new.env(parent = baseenv())
    sys.source(file.path(bindings$root, receipt$stage, "source/tools/profiling/workflow_perf106.R"), helpers)
    pairs <- list()
    valid <- TRUE
    regression <- TRUE
    for (population in unique(vapply(cases, `[[`, character(1), "population"))) for (rep in 1:2) {
        ids <- paste(population, c("control", "candidate"), rep, sep = "-")
        values <- lapply(ids, function(id) {
            x <- receipt$cases[[id]]
            override <- verification_overrides[[id]]
            directory <- if (is.null(override)) file.path(bindings$root, receipt$stage, id) else override$directory
            exit_status <- if (is.null(override)) x$verifier_exit else override$exit_status
            verified <- helpers$readPerf106Verification(directory, dirname(bindings$root), exit_status)
            list(operation = x$operation, verified = verified)
        })
        eligible <- vapply(values, function(x) isTRUE(x$verified$accepted) &&
            identical(x$operation$status, "complete") && isTRUE(x$operation$measurement_valid) &&
            isTRUE(x$operation$lifecycle_passed) && x$operation$peak_bytes < 15.5 * 1024^3 &&
            identical(x$operation$swap_bytes, 0) && x$operation$elapsed_seconds <=
                if (receipt$stage == "runtime-screen") 300 else 1800, logical(1))
        paired <- all(eligible)
        if (paired && population == "p14") paired <- identical(values[[1L]]$verified$value$summary,
            values[[2L]]$verified$value$summary)
        ratios <- if (paired) list(peak = values[[2L]]$operation$peak_bytes / values[[1L]]$operation$peak_bytes,
            elapsed = values[[2L]]$operation$elapsed_seconds / values[[1L]]$operation$elapsed_seconds) else NULL
        within <- paired && ratios$peak <= bindings$thresholds$peak && ratios$elapsed <= bindings$thresholds$elapsed
        pairs[[paste(population, rep, sep = "-")]] <- list(valid = paired, non_regression = within, ratios = ratios)
        valid <- valid && paired
        regression <- regression && within
    }
    list(status = if (!valid) "invalid_or_incomplete" else if (receipt$stage == "runtime-screen")
        "screen_valid" else if (regression) "non_regression_passed" else "non_regression_failed",
        science_outputs_and_measurements_valid = valid, non_regression = regression, pairs = pairs,
        promotion_authorized = FALSE)
}

#' Fingerprint canonical logical pages independently of Parquet row-group layout.
#' @param project,reference Existing immutable output owner and reference.
#' @return Exact page hashes using the original row-name normalization contract.
#' @noRd
fingerprintPerf108LogicalPages <- function(project, reference) {
    identity <- readSavedStoreIdentity(file.path(project, ".workflow-entities/store.json"))
    store <- newWorkflowEntityStore(project, identity$project_id, identity$workflow_id)
    on.exit(closeWorkflowEntityStore(store))
    scope <- newWorkflowTableReadScope(store, reference, reader = "fastduckdb")
    on.exit(releaseWorkflowTableReadScope(scope), add = TRUE, after = FALSE)
    hashes <- character()
    for (offset in seq.int(0L, scope$metadata$rows - 1L, by = 10000L)) {
        data <- readSavedWorkflowTable(store, reference, scope$metadata$columns,
            offset = offset, limit = min(10000L, scope$metadata$rows - offset), scope = scope)$data
        attr(data, "dia_source_order") <- NULL
        attr(data, "row.names") <- .set_row_names(nrow(data))
        hashes[[length(hashes) + 1L]] <- digest::digest(data, algo = "sha256")
    }
    hashes
}
