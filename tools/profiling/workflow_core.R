# CORE series bindings, fixed member schedules, frozen decisions and independent verification.
# Execution stays in the existing PERF coordinator and installed action runner.

#' Frozen per-ticket primary metrics and guards (protocol.md mirrors this table).
#' Every metric is lower-is-better. Metrics compare candidate and control medians;
#' `guards` must hold for each matched pair (control-k vs candidate-k).
#' @noRd
corePrimaryMetrics <- list(
    `006` = list(metrics = "phase:staged_result_read_summary_write", guards = c(peak = 1.00, elapsed = 1.05)),
    `007` = list(metrics = "peak", guards = c(elapsed = 1.05)),
    `008` = list(metrics = c("peak", "phase:relational_validation"), guards = c(elapsed = 1.05)),
    `009` = list(metrics = "peak", guards = c(elapsed = 1.05)),
    `010` = list(metrics = "peak", guards = c(elapsed = 1.05), trade_off = "elapsed"),
    `011` = list(metrics = "peak", guards = c(peak = 1.00, elapsed = 1.05)),
    `012` = list(metrics = character(), guards = numeric(), decision_only = TRUE),
    `013` = list(metrics = c("execution:synchronous_wait_seconds", "execution:sync_forks"),
        guards = c(peak = 1.00, elapsed = 1.05)),
    `014` = list(metrics = c("execution:project_disk_bytes", "execution:head_publishes"),
        guards = c(peak = 1.00, elapsed = 1.05)),
    `017` = list(metrics = "phase_peak:checks", guards = c(peak = 1.00, elapsed = 1.05)))

#' P14 phases every completed CORE member must record (protocol v3).
#' @noRd
coreRequiredPhases <- c("relational_validation", "working_memory_admission",
    "relational_aggregate_order_summary_write", "ordered_result_stage")

#' Worker entry that records phase boundaries beside the request (protocol v3).
#' @param fault `none` or `after_aggregate` (fail once the first ordered result stage has run).
#' @return Rscript expression for the adapters launcher; the scientific call is unchanged.
#' @noRd
corePhaseHookExpression <- function(fault = "none") {
    stopifnot(fault %in% c("none", "after_aggregate"))
    inner <- if (fault == "after_aggregate") paste("function(label, fn) { value <- fn();",
        "if (identical(label, \"ordered_result_stage\"))",
        "stop(\"CORE-004 forced failure after the aggregate stage\"); value }") else "function(label, fn) fn()"
    paste0('local({ d <- dirname(commandArgs(TRUE)[[1L]]); ',
        'events <- file(file.path(d, "phase-events.tsv"), "wt"); ',
        'writeLines("uptime\\tid\\tphase\\tevent", events); flush(events); ',
        'clock <- function() as.numeric(strsplit(readLines("/proc/uptime", n = 1L), " ")[[1L]][[1L]]); ',
        'emit <- function(id, label, event) { writeLines(paste(clock(), id, label, event, sep = "\\t"), events); ',
        'flush(events) }; serial <- 0L; inner <- ', inner, '; ',
        'options(multischolar.dia_rollup_phase = function(label, fn) { serial <<- serial + 1L; id <- serial; ',
        'status <- "failed"; emit(id, label, "start"); on.exit(emit(id, label, status), add = TRUE); ',
        'value <- inner(label, fn); status <- "complete"; value }) }); ',
        'MultiScholaR:::workflowWorkerMain(commandArgs(TRUE))')
}

#' Turn worker phase events and coordinator samples into phase seconds and sampled peaks.
#' @param events Data frame `uptime, id, phase, event` (`start`, then `complete` or `failed`).
#' @param trace Data frame `uptime, current` of the coordinator's 50 ms samples.
#' @return `phase_seconds` (summed intervals), `phase_peak_bytes` (largest sample inside any
#'   interval; NA without a sample), `phase_status` (`failed` > `open` > `complete`).
#' @noRd
corePhaseMeasures <- function(events, trace) {
    stopifnot(is.data.frame(events), all(c("uptime", "id", "phase", "event") %in% names(events)),
        is.data.frame(trace), all(c("uptime", "current") %in% names(trace)))
    starts <- events[events$event == "start", , drop = FALSE]
    ends <- events[events$event %in% c("complete", "failed"), , drop = FALSE]
    labels <- unique(starts$phase)
    seconds <- peaks <- stats::setNames(rep(NA_real_, length(labels)), labels)
    status <- stats::setNames(rep("complete", length(labels)), labels)
    last <- if (nrow(trace)) max(trace$uptime) else max(c(events$uptime, -Inf))
    for (i in seq_len(nrow(starts))) {
        label <- starts$phase[[i]]
        end <- ends[ends$id == starts$id[[i]], , drop = FALSE]
        closed <- nrow(end) == 1L
        finish <- if (closed) end$uptime[[1L]] else last
        state <- if (closed) end$event[[1L]] else "open"
        if (state == "failed" || state == "open" && status[[label]] != "failed") status[[label]] <- state
        seconds[[label]] <- sum(seconds[[label]], finish - starts$uptime[[i]], na.rm = TRUE)
        inside <- trace$current[trace$uptime >= starts$uptime[[i]] & trace$uptime <= finish]
        if (length(inside)) peaks[[label]] <- max(peaks[[label]], inside, na.rm = TRUE)
    }
    list(phase_seconds = as.list(seconds), phase_peak_bytes = as.list(peaks), phase_status = as.list(status))
}

#' Member ceilings in seconds: lifecycle and overflow members are small.
#' @noRd
coreMemberSeconds <- function(stage) if (grepl("^core-(lifecycle-|011-overflow$)", stage)) 300L else 600L

#' Refuse a counted member unless the host meets the plan's start precondition.
#' @param host `list(available, units)` from the coordinator's host reading.
#' @return NULL when ready, otherwise the refusal reason.
#' @noRd
coreSetupRefusal <- function(host) {
    if (!is.numeric(host$available) || length(host$available) != 1L || !is.finite(host$available))
        return("host memory is unobserved")
    if (host$available < 12 * 1024^3) return(sprintf("host MemAvailable %.2f GiB is below 12 GiB",
        host$available / 1024^3))
    if (length(host$units)) return(paste(length(host$units), "other managed unit(s) active"))
    NULL
}

#' Read the tracked CORE bindings into the shape the coordinator and verifier use.
#' @param path `docs/qa/workflow-execution/core/bindings.json`.
#' @return Plain list; no `latest` lookup, every library and input is explicit.
#' @noRd
coreBindings <- function(path) {
    raw <- jsonlite::read_json(path, simplifyVector = TRUE)
    stopifnot(identical(raw$schema, 1L), is.list(raw$harness), is.list(raw$builds))
    input <- function(population) {
        value <- raw$populations[[population]]
        parameters <- value$parameters
        if (!length(parameters)) parameters <- list()
        list(store = value$store, head_sha256 = value$head_sha256, rows = value$rows,
            reference = list(entity_id = value$entity_id), parameters = parameters)
    }
    list(root = raw$harness$root, protocol_sha256 = raw$harness$protocol_sha256,
        builds = raw$builds, outcomes = raw$outcomes, tickets = raw$harness$tickets,
        screen = list(p14 = input("screen_p14")),
        full = list(p14 = input("full_p14"), primary = input("full_primary")),
        p14_oracle = list(screen_pages = raw$oracles$p14$screen_page_hashes,
            full_pages = raw$oracles$p14$full_page_hashes,
            full_summary_json = raw$oracles$p14$full_summary_json,
            receipt = raw$oracles$p14$receipt, receipt_sha256 = raw$oracles$p14$receipt_sha256),
        filter_oracles = list(primary = raw$oracles$filter_primary))
}

#' Select the control build: the most recent ticket whose recorded outcome is accepted.
#' @param outcomes Ordered list of `list(ticket, outcome, build)` records.
#' @return Build key; a reverted or undecided ticket is never a control.
#' @noRd
coreControl <- function(outcomes) {
    accepted <- Filter(\(record) identical(record$outcome, "accepted"), outcomes)
    if (!length(accepted)) stop("no accepted build is recorded")
    accepted[[length(accepted)]]$build
}

#' Expand one CORE stage into its fixed member cases.
#' @param bindings Output of `coreBindings()`.
#' @param stage `core-lifecycle-complete`, `core-lifecycle-fail`, `core-0NN` or `core-011-overflow`.
#' @return Fixed cases in launch order; never an adaptive schedule.
#' @noRd
coreCases <- function(bindings, stage) {
    stopifnot(grepl("^core-(lifecycle-(complete|fail)|0(0[6-9]|1[0-4]|17)|011-overflow)$", stage))
    case <- function(id, role, rep, build_key, input, population, fault = "none", selector) {
        build <- bindings$builds[[build_key]]
        if (is.null(build)) stop("build ", build_key, " is not bound")
        list(id = id, variant = if (population == "p14") "p14" else "physical",
            role = role, control = role == "control", fault = fault, rep = rep,
            population = population, library = build$library, expected_build = build$build,
            exact_sum = build$exact_sum, expect_failure = fault != "none", profile_worker = FALSE,
            phase_events = TRUE,
            binding = digest::digest(list(input, build, stage, rep,
                protocol = bindings$protocol_sha256), algo = "sha256"),
            input = input, selector = selector)
    }
    if (startsWith(stage, "core-lifecycle-")) {
        kind <- sub("core-lifecycle-", "", stage, fixed = TRUE)
        return(list(case("member", "candidate", 1L, "core_004", bindings$screen$p14, "p14",
            fault = if (kind == "fail") "after_aggregate" else "none",
            selector = paste0("core:lifecycle:", kind))))
    }
    if (stage == "core-011-overflow") stop("the overflow member is defined by CORE-011")
    ticket <- sub("core-", "", stage, fixed = TRUE)
    spec <- bindings$tickets[[ticket]]
    if (is.null(spec)) stop("CORE-", ticket, " has no bound control and candidate builds")
    if (ticket == "014") stop("the checkpoint Save/Open vehicle is defined by CORE-014")
    population <- if (ticket == "013") "primary" else "p14"
    input <- bindings$full[[population]]
    roles <- c("control", "candidate", "candidate", "control")
    reps <- c(1L, 1L, 2L, 2L)
    lapply(seq_along(roles), function(i) case(paste(roles[[i]], reps[[i]], sep = "-"), roles[[i]], reps[[i]],
        spec[[roles[[i]]]], input, population,
        selector = paste("core", ticket, roles[[i]], reps[[i]], sep = ":")))
}

#' Map a verifier selector back to its coordinator stage and case id.
#' @noRd
coreSelectorStage <- function(selection) {
    stopifnot(grepl("^core:((00[6-9]|01[0-4]|017):(control|candidate):[12]|lifecycle:(complete|fail)|011:overflow)$", selection))
    parts <- strsplit(sub("^core:", "", selection), ":", fixed = TRUE)[[1L]]
    if (parts[[1L]] == "lifecycle") return(list(stage = paste0("core-lifecycle-", parts[[2L]]), id = "member"))
    if (identical(parts, c("011", "overflow"))) return(list(stage = "core-011-overflow", id = "member"))
    list(stage = paste0("core-", parts[[1L]]), id = paste(parts[[2L]], parts[[3L]], sep = "-"))
}

#' Check one member against the frozen eligibility rules.
#' @param item Coordinator receipt entry with `operation`, `verifier_completion` and `passed`.
#' @param seconds Member ceiling.
#' @return Name of the first failing rule, or NULL when eligible.
#' @noRd
coreIneligible <- function(item, seconds) {
    op <- item$operation
    rules <- list(
        verifier_accepted = isTRUE(item$verifier_completion$accepted),
        status_complete = identical(op$status, "complete"),
        measurement_valid = isTRUE(op$measurement_valid) &&
            !isTRUE(op$execution$monitor_lost_after_result),
        lifecycle_passed = isTRUE(op$lifecycle_passed) && isTRUE(item$passed),
        zero_swap = identical(op$swap_bytes, 0),
        peak_below_early_stop = is.numeric(op$peak_bytes) && length(op$peak_bytes) == 1L &&
            is.finite(op$peak_bytes) && op$peak_bytes > 0 && op$peak_bytes < 15.5 * 1024^3,
        elapsed_within_ceiling = is.numeric(op$elapsed_seconds) && length(op$elapsed_seconds) == 1L &&
            is.finite(op$elapsed_seconds) && op$elapsed_seconds > 0 && op$elapsed_seconds <= seconds)
    failing <- names(rules)[!vapply(rules, isTRUE, logical(1))]
    if (length(failing)) failing[[1L]] else NULL
}

#' Read one named metric from a member's operation receipt.
#' @noRd
coreMetric <- function(op, name) {
    value <- if (name == "peak") op$peak_bytes else if (name == "elapsed") op$elapsed_seconds else
        if (name == "phase_peak:checks") max(unlist(op$phase_peak_bytes[c("relational_validation",
            "working_memory_admission")]) %||% NA_real_) else
        if (startsWith(name, "phase_peak:")) op$phase_peak_bytes[[sub("phase_peak:", "", name, fixed = TRUE)]] else
        if (startsWith(name, "phase:")) op$phase_seconds[[sub("phase:", "", name, fixed = TRUE)]] else
        if (startsWith(name, "execution:")) op$execution[[sub("execution:", "", name, fixed = TRUE)]] else NULL
    if (is.numeric(value) && length(value) == 1L && is.finite(value) && value >= 0) value else NA_real_
}

#' Apply the frozen CORE acceptance rule to one measured ticket's pair of pairs.
#' @param receipt Coordinator stage receipt (`cases` keyed by case id).
#' @param cases Output of `coreCases()` for the stage.
#' @param ticket Three-digit ticket number.
#' @return `accept`, `reject`, `invalid`, `trade_off` or `decision_only`, the first failing rule,
#'   and the median ratios; never promotion or default activation.
#' @noRd
coreDecision <- function(receipt, cases, ticket) {
    spec <- corePrimaryMetrics[[ticket]]
    stopifnot(!is.null(spec))
    result <- list(ticket = ticket, status = "invalid", failing_rule = NULL, ratios = list(),
        pair_ratios = list(), promotion_authorized = FALSE)
    seconds <- coreMemberSeconds(paste0("core-", ticket))
    for (case in cases) {
        rule <- coreIneligible(receipt$cases[[case$id]], seconds)
        if (!is.null(rule)) {
            result$failing_rule <- paste0(case$id, ":", rule)
            return(result)
        }
    }
    if (isTRUE(spec$decision_only)) {
        result$status <- "decision_only"
        return(result)
    }
    median_of <- function(role, name) stats::median(vapply(Filter(\(x) x$role == role, cases),
        \(x) coreMetric(receipt$cases[[x$id]]$operation, name), numeric(1)))
    for (name in spec$metrics) {
        value <- median_of("candidate", name) / median_of("control", name)
        result$ratios[[name]] <- value
        if (is.na(value)) {
            result$failing_rule <- paste0("primary_metric_unavailable:", name)
            return(result)
        }
    }
    # Guards hold per matched pair (protocol v2, user 2026-10-08): no pair may breach a limit.
    guards <- spec$guards
    pair_ratio <- function(name, rep) coreMetric(receipt$cases[[paste0("candidate-", rep)]]$operation, name) /
        coreMetric(receipt$cases[[paste0("control-", rep)]]$operation, name)
    reps <- sort(unique(vapply(cases, `[[`, integer(1), "rep")))
    for (name in names(guards)) result$pair_ratios[[name]] <- vapply(reps, \(rep) pair_ratio(name, rep), numeric(1))
    worse <- spec$metrics[vapply(spec$metrics, \(name) result$ratios[[name]] > 1, logical(1))]
    breached <- names(guards)[vapply(names(guards), \(name) any(!(result$pair_ratios[[name]] <= guards[[name]])),
        logical(1))]
    result$status <- "reject"
    if (length(worse)) {
        result$failing_rule <- paste0("primary_metric_worse:", worse[[1L]])
    } else if (length(breached)) {
        result$failing_rule <- paste0("guard_breached:", breached[[1L]])
        if (identical(breached, spec$trade_off)) result$status <- "trade_off"
    } else {
        result$status <- "accept"
    }
    result
}

#' Classify a lifecycle member; a failure member can never be relabelled as success.
#' @param receipt,cases Coordinator stage receipt and its single fixed case.
#' @return `lifecycle_complete_observed`, `lifecycle_failure_observed` or `invalid` with the failing rule.
#' @noRd
coreLifecycleDecision <- function(receipt, cases) {
    stopifnot(length(cases) == 1L)
    case <- cases[[1L]]
    item <- receipt$cases[[case$id]]
    op <- item$operation
    result <- list(status = "invalid", failing_rule = NULL, promotion_authorized = FALSE,
        expect_failure = case$expect_failure)
    status <- op$phase_status
    if (!case$expect_failure) {
        result$failing_rule <- coreIneligible(item, coreMemberSeconds("core-lifecycle-complete"))
        # Protocol v3: phase timing must observe the work through completion.
        if (is.null(result$failing_rule) && !(all(coreRequiredPhases %in% names(status)) &&
            all(unlist(status) == "complete") && all(is.finite(unlist(op$phase_peak_bytes[coreRequiredPhases])))))
            result$failing_rule <- "phase_timing"
        if (is.null(result$failing_rule)) result$status <- "lifecycle_complete_observed"
        return(result)
    }
    rules <- list(verifier_accepted = isTRUE(item$verifier_completion$accepted),
        failure_recorded = identical(op$status, "worker_failed"),
        lifecycle_passed = isTRUE(op$lifecycle_passed) && isTRUE(item$passed),
        zero_swap = identical(op$swap_bytes, 0),
        cleanup = isTRUE(op$execution$cleanup$worker_stopped) && isTRUE(op$execution$cleanup$monitor_stopped),
        phase_timing = identical(status$ordered_result_stage, "failed") &&
            identical(status$relational_validation, "complete"))
    failing <- names(rules)[!vapply(rules, isTRUE, logical(1))]
    if (length(failing)) result$failing_rule <- failing[[1L]] else result$status <- "lifecycle_failure_observed"
    result
}

#' Verify one CORE member outside its measured application slice.
#' @param run_dir,repo Guarded verifier output and frozen source snapshot.
#' @param selection `core:` verifier selector naming the member.
#' @return Source-bound verification receipt; no measured work is repeated.
#' @noRd
verifyCoreMember <- function(run_dir, repo, selection) {
    target <- coreSelectorStage(selection)
    bindings <- coreBindings(file.path(repo, "docs/qa/workflow-execution/core/bindings.json"))
    cases <- coreCases(bindings, target$stage)
    case <- Filter(\(x) x$id == target$id, cases)[[1L]]
    owned <- file.path(bindings$root, target$stage, case$id)
    path <- file.path(owned, "verification.rds")
    stopifnot(!file.exists(path))
    receipt <- list(status = "starting", run_dir = run_dir, member_id = case$id, selector = selection)
    on.exit({
        if (identical(receipt$status, "starting")) receipt$status <- "failed"
        saveRDS(receipt, path)
    }, add = TRUE)
    operation_path <- file.path(owned, "operation-receipt.rds")
    receipt$operation_sha256 <- digest::digest(file = operation_path, algo = "sha256")
    operation <- readRDS(operation_path)
    member <- readRDS(file.path(owned, "member-result.rds"))
    expected_status <- if (case$expect_failure) "worker_failed" else "complete"
    stopifnot(identical(operation$status, expected_status), identical(member$status, expected_status),
        identical(operation$config$build, case$expected_build), identical(operation$binding, case$binding),
        identical(member$project_root, file.path(owned, "project")),
        identical(operation$member$input_entity_id, case$input$reference$entity_id),
        identical(digest::digest(file = file.path(case$input$store, "dia-state-head.rds"),
            algo = "sha256"), case$input$head_sha256))
    project <- member$project_root
    owner <- readSavedStoreIdentity(file.path(project, ".workflow-entities/store.json"))
    manager <- newDIAStateAdapter(newWorkflowEntityStore(project, owner$project_id, owner$workflow_id), reopen = TRUE)
    on.exit(closeDIAStateAdapter(manager), add = TRUE, after = FALSE)
    checkpoint <- readDIAStateMetadata(manager)$data
    request <- file.path(project, ".workflow-entities/.workflow-export", operation$execution$request_id)
    journal <- jsonlite::read_json(file.path(request, "journal.json"))
    completion <- readRDS(file.path(request, "completion.rds"))
    worker <- readRDS(file.path(request, "worker-result.rds"))
    # CORE-003 states: a conflict, an invalid guard payload or a lost monitor never verifies.
    stopifnot(!identical(journal$status, "conflict"), !identical(completion$status, "guard_payload_invalid"),
        !isTRUE(completion$monitor_lost_after_result), identical(completion$status, expected_status),
        isTRUE(completion$cleanup$worker_stopped), isTRUE(completion$cleanup$monitor_stopped))
    if (case$expect_failure) {
        stopifnot(identical(journal$status, "worker_failed"), identical(journal$completion_status, "worker_failed"),
            identical(worker$status, "failed"),
            grepl("CORE-004 forced failure after the aggregate stage", worker$error, fixed = TRUE),
            identical(checkpoint$components$peptide_data$entity_id, case$input$reference$entity_id))
        receipt$failure <- list(worker_error = worker$error, journal_status = journal$status,
            head_entity = checkpoint$components$peptide_data$entity_id, completion_status = completion$status)
        receipt$status <- "complete"
        return(receipt)
    }
    stopifnot(identical(journal$status, "complete"), identical(case$population, "p14"))
    execution <- worker$result$execution
    screen <- identical(case$input$rows, bindings$screen$p14$rows)
    reference <- checkpoint$components$peptide_data
    hashes <- fingerprintPerf108LogicalPages(project, reference)
    expected <- bindings$p14_oracle[[if (screen) "screen_pages" else "full_pages"]]
    encode <- function(value) jsonlite::toJSON(value, auto_unbox = TRUE, null = "null", digits = NA)
    summary <- as.character(encode(checkpoint$slots$args$rollUpPrecursorToPeptide$rollup_summary))
    figures <- list.files(file.path(request, "required-output"),
        "4_precursor_rollup_.*[.](png|pdf)$", recursive = TRUE, full.names = TRUE)
    figures <- figures[!grepl("/figures/", figures)]
    stopifnot(identical(unname(hashes), unname(expected)),
        reference$metadata$shape[[1L]] == if (screen) 480000L else 12000000L,
        identical(execution$exact_sum_sha256, case$exact_sum),
        identical(execution$query_owner, "fastduckdb_c"),
        identical(execution$profile$name, "P14"), identical(execution$engine_threads, 8L),
        identical(execution$writer_backend, "fast_r"), identical(execution$writer_threads, 8L),
        !is.null(execution$admission$worker_cap), length(figures) == 12L,
        all(file.size(figures) > 0), is.null(assertDIAExportReady(manager)))
    if (!screen) stopifnot(identical(summary, bindings$p14_oracle$full_summary_json))
    receipt$summary <- summary
    receipt$pages <- hashes
    receipt$pages_match <- TRUE
    receipt$figures <- stats::setNames(file.size(figures), basename(figures))
    receipt$execution <- execution
    receipt$status <- "complete"
    receipt
}
