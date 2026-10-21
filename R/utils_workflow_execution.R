# Internal WFX execution boundary. Scientific adapters are deliberately finite.

#' Evaluate frozen paired qualification evidence without executing an operation.
#' @param evidence Bounded typed qualification record (schema 1).
#' @param expected_bindings Independently frozen operation/input/runtime/boundary identities.
#' @param limits Server-owned absolute caps; NULL uses the common qualification maximum.
#' @return One deterministic decision; malformed or incomplete evidence is held.
#' @noRd
evaluateWorkflowQualification <- function(evidence, expected_bindings, limits = NULL) {
    if (is.list(evidence) && identical(evidence$schema, 2L))
        return(evaluatePerfWorkflowQualification(evidence, expected_bindings, limits))
    held <- list(status = "held", ratios = NULL)
    flags <- c("domain_supported", "science", "lifecycle", "integrity", "outputs",
        "absolute_limits", "measurement_valid", "bindings_match", "runtime_eligible",
        "predecessor_benefit_bound")
    binding_names <- c("operation", "input", "parameters", "runtime", "boundary", "policy", "outputs")
    valid_bindings <- function(x) is.list(x) && identical(sort(names(x)), sort(binding_names)) &&
        all(vapply(x, function(value) is.character(value) && length(value) == 1L &&
            !is.na(value) && grepl("^[a-f0-9]{64}$", value), logical(1)))
    if (!is.list(evidence) || !identical(evidence$schema, 1L) ||
        !valid_bindings(expected_bindings) || !valid_bindings(evidence$bindings) ||
        !identical(evidence$bindings, expected_bindings) ||
        !all(vapply(flags, function(name) is.logical(evidence[[name]]) &&
            length(evidence[[name]]) == 1L && !is.na(evidence[[name]]), logical(1))) ||
        !is.list(evidence$pairs) || length(evidence$pairs) != 2L ||
        !is.character(evidence$open_holds) || anyNA(evidence$open_holds)) return(held)
    status <- function(value) list(status = value, ratios = NULL)
    for (pair in evidence$pairs) {
        if (!is.list(pair)) return(held)
        for (arm in c("control", "candidate")) {
            member <- pair[[arm]]
            if (!is.list(member) || !is.character(member$status) || length(member$status) != 1L ||
                is.na(member$status)) return(held)
        }
    }
    if (!evidence$domain_supported) return(status("unsupported"))
    if (!all(unlist(evidence[c("science", "lifecycle", "integrity", "outputs", "absolute_limits")])))
        return(status("rejected"))
    known <- c("complete", "memory_guard", "time_guard", "dependency", "science", "publication")
    if (any(vapply(evidence$pairs, function(pair)
        !pair$control$status %in% known || !pair$candidate$status %in% known, logical(1)))) return(held)
    if (any(vapply(evidence$pairs, function(pair) pair$candidate$status != "complete" ||
        pair$control$status %in% c("science", "publication"), logical(1)))) return(status("rejected"))
    if (!is.character(evidence$goal) || length(evidence$goal) != 1L || is.na(evidence$goal) ||
        !evidence$goal %in% c("benefit", "integration_non_regression", "capacity_only") ||
        length(evidence$open_holds) || !evidence$measurement_valid || !evidence$bindings_match ||
        !evidence$runtime_eligible || evidence$goal == "integration_non_regression" &&
        !evidence$predecessor_benefit_bound) return(held)
    positive <- function(value) is.numeric(value) && length(value) == 1L &&
        !is.na(value) && is.finite(value) && value > 0
    counters <- function(member) positive(member$peak_bytes) && positive(member$elapsed_seconds)
    if (is.null(limits)) limits <- list(memory_max = 16 * 1024^3, deadline = 1800)
    if (!is.list(limits) || !positive(limits$memory_max) || !positive(limits$deadline) ||
        limits$memory_max > 16 * 1024^3 || limits$deadline > 1800) return(held)
    for (pair in evidence$pairs) {
        if (counters(pair$candidate) && (pair$candidate$peak_bytes > limits$memory_max ||
            pair$candidate$elapsed_seconds > limits$deadline)) return(status("rejected"))
        if (pair$control$status == "complete" && counters(pair$control) &&
            (pair$control$peak_bytes > limits$memory_max || pair$control$elapsed_seconds > limits$deadline)) return(held)
    }
    if (evidence$goal == "capacity_only") {
        good <- vapply(evidence$pairs, function(pair) pair$control$status %in%
            c("memory_guard", "time_guard") && isTRUE(pair$control$input_valid) &&
            counters(pair$candidate), logical(1))
        return(if (all(good)) status("capacity_only") else held)
    }
    if (!all(vapply(evidence$pairs, function(pair) pair$control$status == "complete" &&
        counters(pair$control) && counters(pair$candidate), logical(1)))) return(held)
    ratios <- lapply(evidence$pairs, function(pair) list(
        peak = pair$candidate$peak_bytes / pair$control$peak_bytes,
        elapsed = pair$candidate$elapsed_seconds / pair$control$elapsed_seconds))
    limit <- if (evidence$goal == "benefit") 0.90 else 1.05
    passed <- vapply(ratios, function(pair) pair$peak <= limit && pair$elapsed <= 1.20, logical(1))
    list(status = if (all(passed)) evidence$goal else if (any(passed)) "inconclusive" else "rejected",
        ratios = ratios)
}

#' Read a bounded internal execution journal using existing publication paths.
#' @param paths Owned workflowExportPaths result.
#' @return Journal or NULL; corruption is never interpreted as an empty journal.
#' @noRd
readWorkflowExecutionJournal <- function(paths) {
    path <- file.path(paths$request_dir, "journal.json")
    if (!file.exists(path)) return(NULL)
    if (file.size(path) > 1048576 || nzchar(Sys.readlink(path)))
        cli::cli_abort("execution journal is invalid", class = "workflow_execution_invalid")
    journal <- jsonlite::read_json(path)
    if (!identical(journal$schema, "multischolar.workflow_action.v1"))
        cli::cli_abort("execution journal schema differs", class = "workflow_execution_invalid")
    journal
}

#' Identify installed execution code; development/source injection cannot qualify it.
#' @return Immutable installed application binary identity.
#' @noRd
workflowExecutionBuildIdentity <- function() {
    path <- system.file("R", "MultiScholaR.rdb", package = "MultiScholaR")
    if (!nzchar(path) || !file.exists(path)) cli::cli_abort(
        "installed execution runtime required", class = "workflow_execution_unsupported")
    digest::digest(file = path, algo = "sha256")
}

#' Reject stateful work inside a worker before creating state artifacts.
#' @return NULL invisibly for a controller/legacy caller.
#' @noRd
assertWorkflowController <- function() {
    if (identical(Sys.getenv("MULTISCHOLAR_WORKFLOW_WORKER"), "1"))
        cli::cli_abort("workers cannot mutate workflow heads", class = "workflow_execution_worker_mutation")
    invisible(NULL)
}

#' Acquire the project mutation lease, including for legacy synchronous publication.
#' @param store Existing entity owner.
#' @return A temporary OS lock or NULL when this controller already holds it.
#' @noRd
acquireWorkflowMutationLease <- function(store) {
    assertWorkflowController()
    checkWorkflowEntityStore(store)
    owner <- store$execution_owner
    if (!is.null(owner)) {
        if (owner$closed || !identical(owner$pid, Sys.getpid()) || is.null(owner$lease))
            cli::cli_abort("execution controller lease is unavailable", class = "workflow_execution_busy")
        return(NULL)
    }
    control_path <- file.path(store$directory, ".workflow-execution-control.rds")
    if (file.exists(control_path)) {
        if (file.size(control_path) > 65536) cli::cli_abort("execution control metadata exceeds limit")
        control <- readRDS(control_path)
        if (isTRUE(control$active) && identical(control$incarnation, store$identity$incarnation) &&
            identical(control$process, workflowExecutionProcess()))
            cli::cli_abort("this project already has a managed controller", class = "workflow_execution_busy")
    }
    lock <- filelock::lock(file.path(store$directory, ".workflow-mutation.lock"), timeout = 0)
    if (is.null(lock)) cli::cli_abort("another controller owns this project", class = "workflow_execution_busy")
    lock
}

#' Persist only bounded coordination metadata in the existing store.
#' @param store Existing owner.
#' @param control Internal lease/generation record, never scientific data.
#' @return NULL invisibly.
#' @noRd
writeWorkflowExecutionControl <- function(store, control) {
    checkWorkflowEntityStore(store)
    if (as.numeric(object.size(control)) > 65536) cli::cli_abort("oversized execution control")
    path <- file.path(store$directory, ".workflow-execution-control.rds")
    temp <- tempfile("execution-control-", tmpdir = store$directory)
    on.exit(if (file.exists(temp)) unlink(temp))
    saveRDS(control, temp, compress = FALSE, version = 3)
    syncSavedWorkflowPaths(store, c(temp, store$directory))
    if (!file.rename(temp, path)) cli::cli_abort("execution control publication failed")
    syncSavedWorkflowPaths(store, c(path, store$directory))
    invisible(NULL)
}

#' Observe a durable head by immutable identity without hydrating its state.
#' @param store Existing project owner.
#' @return Descriptor or NULL before initial publication.
#' @noRd
workflowExecutionHead <- function(store) {
    path <- file.path(store$directory, "dia-state-head.rds")
    if (!file.exists(path)) return(NULL)
    if (file.size(path) > 65536) cli::cli_abort("oversized workflow head")
    readRDS(path)
}

#' Prepare journal reconciliation before the existing atomic head advancement.
#' @param context Existing DIA owner.
#' @param reference Already durable proposed head.
#' @return NULL invisibly.
#' @noRd
prepareWorkflowHeadAcceptance <- function(context, reference) {
    owner <- context$store$execution_owner
    if (is.null(owner) || is.null(owner$accepting)) return(invisible(NULL))
    paths <- workflowExportPaths(context$store$directory, owner$accepting)
    journal <- readWorkflowExecutionJournal(paths)
    pending <- readRDS(file.path(paths$request_dir, "pending-result.rds"))
    snapshot <- hydrateSavedWorkflowCheckpoint(context$store, reference)
    accepted <- list(reference = pending$reference, summary = pending$summary,
        dependencies = pending$dependencies, state_reference = snapshot$states[[snapshot$current_state]],
        operation = pending$operation %||% workflowFixtureOperation)
    result_path <- file.path(paths$request_dir, "accepted-result.rds")
    writeWorkflowExecutionMessage(result_path, accepted)
    syncSavedWorkflowPaths(context$store, result_path)
    journal$result_sha256 <- digest::digest(file = result_path, algo = "sha256")
    journal$status <- "committing"
    journal$proposed_head <- reference$entity_id
    journal$proposed_head_path <- reference$payload$relative_path
    # Typed descriptor is separately durable; JSON is never substituted for its R types.
    proposed <- file.path(paths$request_dir, "proposed-head.rds")
    saveRDS(reference, proposed, compress = FALSE, version = 3)
    syncSavedWorkflowPaths(context$store, c(proposed, paths$request_dir))
    journal$proposed_descriptor_sha256 <- digest::digest(file = proposed, algo = "sha256")
    writeWorkflowExportJournal(paths, journal)
    invisible(NULL)
}

#' Record a completed head update and invalidate earlier dispatch generations.
#' @param context Existing DIA owner.
#' @return NULL invisibly.
#' @noRd
finishWorkflowHeadAcceptance <- function(context) {
    owner <- context$store$execution_owner
    if (is.null(owner)) return(invisible(NULL))
    owner$control$generation <- owner$control$generation + 1L
    writeWorkflowExecutionControl(context$store, owner$control)
    if (!is.null(owner$accepting)) {
        paths <- workflowExportPaths(context$store$directory, owner$accepting)
        journal <- readWorkflowExecutionJournal(paths)
        journal$status <- "science_committed"
        journal$accepted_head <- context$head$entity_id
        writeWorkflowExportJournal(paths, journal)
    }
    invisible(NULL)
}

#' Process-owned controller for one existing project store.
#' @noRd
WorkflowExecutionController <- R6::R6Class("WorkflowExecutionController", cloneable = FALSE,
    public = list(manager = NULL, application = NULL, lease = NULL, control = NULL,
        pid = NULL, closed = FALSE, accepting = NULL, job = NULL))

#' Acquire a project controller and reconcile an interrupted atomic head publication.
#' @param manager Existing reference-state adapter.
#' @param application Validated managed application owner.
#' @return Controller holding the OS mutation lease until explicit close.
#' @noRd
newWorkflowExecutionController <- function(manager, application) {
    assertWorkflowController()
    checkWorkflowApplication(application)
    store <- manager$.context$store
    lease <- acquireWorkflowMutationLease(store)
    if (is.null(lease)) cli::cli_abort("controller already registered", class = "workflow_execution_busy")
    ready <- FALSE
    on.exit(if (!ready) filelock::unlock(lease))
    path <- file.path(store$directory, ".workflow-execution-control.rds")
    control <- if (file.exists(path)) readRDS(path) else list(schema = 1L,
        incarnation = store$identity$incarnation, generation = 0L, active_request = NULL)
    if (!identical(control$schema, 1L) || !identical(control$incarnation, store$identity$incarnation) ||
        !is.numeric(control$generation) || length(control$generation) != 1L ||
        !is.finite(control$generation) || control$generation < 0)
        cli::cli_abort("execution control differs", class = "workflow_execution_invalid")
    owner <- WorkflowExecutionController$new()
    owner$manager <- manager
    owner$application <- application
    owner$lease <- lease
    owner$pid <- Sys.getpid()
    owner$control <- control
    recoverWorkflowHeadAcceptance(owner)
    owner$control$generation <- owner$control$generation + 1L
    owner$control$pid <- Sys.getpid()
    owner$control$process <- workflowExecutionProcess()
    owner$control$active <- TRUE
    owner$control$controller_id <- digest::digest(list(Sys.getpid(), Sys.time(), tempfile()), algo = "sha256")
    writeWorkflowExecutionControl(store, owner$control)
    store$execution_owner <- owner
    application$controllers[[store$directory]] <- owner
    ready <- TRUE
    owner
}

#' Reconcile one journaled head intent before admitting new mutations.
#' @param owner New leased controller; no workers have been admitted yet.
#' @return NULL invisibly; unknown heads hold recovery rather than guessing.
#' @noRd
recoverWorkflowHeadAcceptance <- function(owner) {
    request <- owner$control$active_request
    if (is.null(request)) return(invisible(NULL))
    context <- owner$manager$.context
    paths <- workflowExportPaths(context$store$directory, request)
    journal <- readWorkflowExecutionJournal(paths)
    if (is.null(journal)) cli::cli_abort("active execution journal is missing")
    refuse <- function() cli::cli_abort(c("head recovery conflicts with a later mutation",
        i = "Verify the durable head, then call {.code MultiScholaR:::acknowledgeWorkflowExecutionConflict(root, \"{request}\", head)}."),
        class = c("workflow_execution_conflict", "workflow_execution_stale"), request_id = request)
    if (identical(journal$status, "conflict")) refuse()
    if (identical(journal$status, "committing")) {
        descriptor <- file.path(paths$request_dir, "proposed-head.rds")
        if (!file.exists(descriptor) || file.size(descriptor) > 65536 ||
            !identical(digest::digest(file = descriptor, algo = "sha256"), journal$proposed_descriptor_sha256))
            cli::cli_abort("proposed head is invalid", class = "workflow_execution_invalid")
        proposed <- readRDS(descriptor)
        validateSavedWorkflowReference(context$store, proposed, payload = TRUE)
        current <- workflowExecutionHead(context$store)
        if (identical(current$entity_id, journal$expected_head)) {
            temporary <- tempfile("execution-recover-", tmpdir = context$store$directory)
            saveRDS(proposed, temporary, compress = FALSE, version = 3)
            syncSavedWorkflowPaths(context$store, temporary)
            if (!file.rename(temporary, context$head_path)) cli::cli_abort("head recovery failed")
            syncSavedWorkflowPaths(context$store, c(context$head_path, context$store$directory))
            current <- proposed
        }
        if (!identical(current$entity_id, proposed$entity_id)) {
            journal$status <- "conflict"
            journal$conflict_head <- current$entity_id
            writeWorkflowExportJournal(paths, journal)
            refuse()
        }
        journal$status <- "science_committed"
        journal$accepted_head <- proposed$entity_id
        writeWorkflowExportJournal(paths, journal)
    }
    current <- workflowExecutionHead(context$store)
    if (!is.null(current)) {
        context$head <- current
        context$snapshot <- resolveDIAStateValue(context, current)
    }
    if (!journal$status %in% c("science_committed", "publication_pending", "complete")) {
        journal$status <- "controller_lost"
        writeWorkflowExportJournal(paths, journal)
    }
    owner$control$active_request <- NULL
    invisible(NULL)
}

#' Release a recovery conflict after an operator has verified the durable head.
#' @param root Project directory containing `.workflow-entities`.
#' @param request_id Request recorded as `conflict` in its journal.
#' @param head Entity id the operator verified as the current durable head, or NULL
#'   when the project has no durable head.
#' @return The acknowledged journal, invisibly; the proposed head is never applied.
#' @noRd
acknowledgeWorkflowExecutionConflict <- function(root, request_id, head) {
    identity <- readSavedStoreIdentity(file.path(root, ".workflow-entities", "store.json"))
    store <- newWorkflowEntityStore(root, identity$project_id, identity$workflow_id)
    on.exit(closeWorkflowEntityStore(store), add = TRUE)
    lease <- acquireWorkflowMutationLease(store)
    if (is.null(lease)) cli::cli_abort("controller already registered", class = "workflow_execution_busy")
    on.exit(filelock::unlock(lease), add = TRUE, after = FALSE)
    control <- readRDS(file.path(store$directory, ".workflow-execution-control.rds"))
    paths <- workflowExportPaths(store$directory, request_id)
    journal <- readWorkflowExecutionJournal(paths)
    if (!identical(control$active_request, request_id) || !identical(journal$status, "conflict"))
        cli::cli_abort("request {request_id} holds no recovery conflict", class = "workflow_execution_invalid")
    valid_head <- is.null(head) || is.character(head) && length(head) == 1L && !is.na(head)
    if (!valid_head || !identical(workflowExecutionHead(store)$entity_id, head))
        cli::cli_abort("acknowledged head is not the durable head", class = "workflow_execution_stale")
    journal$acknowledged_head <- head
    writeWorkflowExportJournal(paths, journal)
    control$active_request <- NULL
    writeWorkflowExecutionControl(store, control)
    invisible(journal)
}

#' Cancel before acknowledging Revert/project replacement or session loss.
#' @param owner Project execution owner.
#' @param reason Bounded cancellation classification.
#' @return NULL invisibly after worker termination request and generation advance.
#' @noRd
cancelWorkflowAction <- function(owner, reason = "cancelled") {
    if (owner$closed) return(invisible(NULL))
    owner$control$generation <- owner$control$generation + 1L
    if (is.null(owner$job)) owner$control$active_request <- NULL
    if (workflowExecutionStoreValid(owner$manager$.context$store))
        writeWorkflowExecutionControl(owner$manager$.context$store, owner$control)
    if (!is.null(owner$job)) {
        job <- owner$job
        job$cancelled <- reason
        if (!is.null(job$worker_unit)) workflowSystemd(c("stop", job$worker_unit), required = FALSE)
        result <- finishWorkflowExecutionJob(owner, reason)
        if (identical(result$status, "cleanup_failed")) cli::cli_abort(
            "worker cleanup is incomplete; controller ownership retained",
            class = "workflow_execution_cleanup_failed")
    }
    invisible(NULL)
}

#' Release only the controller lease and owned execution resources.
#' @param owner Project controller.
#' @return NULL invisibly; durable artifacts/journals remain readable.
#' @noRd
closeWorkflowExecutionController <- function(owner) {
    if (owner$closed) return(invisible(NULL))
    store <- owner$manager$.context$store
    if (!is.null(owner$job)) cancelWorkflowAction(owner, "controller_closed")
    owner$control$active <- FALSE
    owner$control$active_request <- NULL
    if (workflowExecutionStoreValid(store)) writeWorkflowExecutionControl(store, owner$control)
    filelock::unlock(owner$lease)
    owner$lease <- NULL
    store$execution_owner <- NULL
    owner$application$controllers[[store$directory]] <- NULL
    owner$closed <- TRUE
    invisible(NULL)
}

#' Run a bounded systemd command with argument-vector quoting, never shell interpolation.
#' @param args Fixed command arguments and validated internal identifiers.
#' @param required Whether command failure is an unsupported-runtime error.
#' @param command System command, selected only by package code.
#' @return Command output with its exit status.
#' @noRd
workflowSystemd <- function(args, required = TRUE, command = "systemctl") {
    binary <- Sys.which(command)
    if (!nzchar(binary)) cli::cli_abort("managed execution command unavailable",
        class = "workflow_execution_unsupported")
    output <- suppressWarnings(system2(binary, shQuote(c("--user", args)), stdout = TRUE, stderr = TRUE))
    status <- attr(output, "status") %||% 0L
    if (required && status != 0L) cli::cli_abort("managed execution command failed: {paste(output, collapse = ' ')}",
        class = "workflow_execution_unsupported")
    attr(output, "exit_status") <- status
    output
}

#' Read a fixed set of cached service properties, including after cgroup teardown.
#' @param unit Internally generated service/slice identifier.
#' @return Named character properties; missing services remain explicit.
#' @noRd
workflowServiceProperties <- function(unit) {
    if (!is.character(unit) || length(unit) != 1L ||
        !grepl("^multischolarwfx[a-f0-9]+(-[a-z0-9]+)?[.](service|slice)$", unit))
        cli::cli_abort("invalid managed unit", class = "workflow_execution_invalid")
    output <- workflowSystemd(c("show", unit,
        "--property=ActiveState,SubState,Result,MainPID,ControlGroup,MemoryPeak,MemorySwapPeak,CPUUsageNSec,IOReadBytes,IOWriteBytes,MemoryMax,MemorySwapMax"), required = FALSE)
    if (attr(output, "exit_status") != 0L) return(list(available = FALSE))
    fields <- strsplit(output[grepl("=", output, fixed = TRUE)], "=", fixed = TRUE)
    out <- stats::setNames(lapply(fields, function(x) paste(x[-1L], collapse = "=")),
        vapply(fields, `[[`, character(1), 1L))
    out$available <- TRUE
    out
}

#' Read one kernel memory counter without guessing missing values.
#' @param path Counter file.
#' @return Nonnegative finite bytes.
#' @noRd
workflowMemoryCounter <- function(path) {
    value <- tryCatch(readLines(path, n = 1L, warn = FALSE), error = function(error) character())
    if (length(value) != 1L || !grepl("^[0-9]+$", value)) cli::cli_abort(
        "memory observation is unavailable", class = "workflow_execution_unsupported")
    as.numeric(value)
}

#' Read available host memory from the kernel, including reclaimable cache.
#' @return Available bytes; unsupported/missing counters fail closed.
#' @noRd
workflowHostAvailable <- function() {
    line <- grep("^MemAvailable:", readLines("/proc/meminfo", warn = FALSE), value = TRUE)
    fields <- if (length(line) == 1L) strsplit(trimws(line), "[[:space:]]+")[[1L]] else character()
    value <- if (length(fields) >= 2L) suppressWarnings(as.numeric(fields[[2L]])) * 1024 else NA_real_
    if (length(value) != 1L || !is.finite(value) || value < 0) cli::cli_abort(
        "host memory unavailable: /proc/meminfo has no usable MemAvailable (Linux 3.14 or later required)",
        class = "workflow_execution_unsupported")
    value
}

#' Return the server-owned policy; small diagnostics never qualify scientific inputs.
#' @param diagnostic Explicit internal lifecycle-fixture mode.
#' @return Fixed resource envelope, never populated from GUI input.
#' @noRd
workflowExecutionPolicy <- function(diagnostic = FALSE) {
    if (diagnostic) return(list(schema = 1L, diagnostic = TRUE, threads = 8L,
        memory_max = 4 * 1024^3, early_stop = 2.9 * 1024^3, swap_max = 0,
        growth_reserve = 128 * 1024^2, startup_headroom = 2 * 1024^3,
        host_floor = 2 * 1024^3, engine_limit = 12 * 1024^3, deadline = 90L,
        monitor_interval = 0.05))
    list(schema = 1L, diagnostic = FALSE, threads = 8L, memory_max = 16 * 1024^3,
        early_stop = 15.5 * 1024^3, swap_max = 0, growth_reserve = 1024^3,
        startup_headroom = 4 * 1024^3, host_floor = 2 * 1024^3,
        engine_limit = 12 * 1024^3, deadline = 1800L, monitor_interval = 0.05)
}

#' Bound the worker by the slice and by what the host has free at dispatch (CORE-010).
#' @param policy Server-owned execution policy.
#' @param current Slice memory already in use (controller and monitor).
#' @param host Host available memory at dispatch.
#' @return `cap` in bytes, its `source` (`slice` or `host`) and the slice `reservation`.
#' @noRd
workflowAdmissionCap <- function(policy, current, host) {
    if (!is.numeric(current) || length(current) != 1L || !is.finite(current) ||
        !is.numeric(host) || length(host) != 1L || !is.finite(host))
        cli::cli_abort("admission memory readings unavailable", class = "workflow_execution_unsupported")
    reservation <- current + policy$growth_reserve
    slice <- policy$memory_max - reservation
    # The monitor stops the action below the host floor, so the growth reserve also stays free
    # above it (CORE-S3-P03); other programs keep what they hold.
    free <- host - policy$host_floor - policy$growth_reserve
    list(cap = min(slice, free), source = if (free < slice) "host" else "slice", reservation = reservation)
}

#' Managed application resource owner; no scientific state is stored here.
#' @noRd
WorkflowExecutionApplication <- R6::R6Class("WorkflowExecutionApplication", cloneable = FALSE,
    public = list(config = NULL, cgroup = NULL, controllers = list(), active = NULL,
        sessions = 0L, pid = NULL, closed = FALSE, entry_hook = NULL))

#' The dispatch rule (KIT-008, audit B18): a module runs a registered operation on the managed worker
#' when its Shiny session carries a `WorkflowExecutionApplication` (set by `workflowManagedServer`
#' under `launchWorkflowApplication`: Linux, systemd); otherwise the same declaration runs in process.
#' @param session Current Shiny session.
#' @return The managed application, or NULL for an unmanaged launch (`run_app()`, macOS, Windows).
#' @noRd
workflowSessionApplication <- function(session) {
    application <- if (is.null(session)) NULL else session$userData$workflow_execution
    if (inherits(application, "WorkflowExecutionApplication") && !isTRUE(application$closed)) application else NULL
}

#' The completion context a GUI step hands to the managed action (KIT-008, audit PE-02): the bus's
#' accumulated peptide QC parameters and the prior filtering progress, so the controller's completion
#' writes the same `qc_progress` product the in-process path writes.
#' @param workflow_data Workflow bus.
#' @return List `peptide_qc` (list) and `prior_progress` (slot list or NULL).
#' @noRd
workflowCompletionContext <- function(workflow_data) {
    qc_params <- shiny::isolate(workflow_data$qc_params)
    prior <- shiny::isolate(workflow_data$filtering_progress)
    list(peptide_qc = qc_params$peptide_qc %||% list(),
        prior_progress = if (methods::is(prior, "FilteringProgress")) stats::setNames(
            lapply(methods::slotNames(prior), function(name) methods::slot(prior, name)), methods::slotNames(prior)) else NULL)
}

#' The completion context a managed step hands over: the peptide QC completion role reads it (the bus's
#' peptide QC parameters and the prior filtering progress); every other role gets none. KIT-017 counted pair:
#' the DA publication step handed over the 300-sample filtering progress, which exceeds the 4 MiB bound, and the
#' step failed at full size while the screen corpus passed.
#' @noRd
workflowManagedCompletionContext <- function(operation, workflow_data) {
    entry <- workflowOperationEntry(operation)
    if (identical(entry$completion, "peptide")) workflowCompletionContext(workflow_data) else NULL
}

#' The browser's ready message for a managed step: the accepted text and, when the step shows a plot at
#' completion, the plot's output id. A step whose completion shows no plot (`ui$plot` NULL: the DA publication,
#' whose static volcano renders only when its sub-tab is shown; KIT-017 audit PH-04, KIT-037) is acknowledged on
#' its text alone and the browser reports the absent plot as ready.
#' @noRd
workflowManagedReadyPayload <- function(ns, action, ui) {
    list(request_id = action$reservation$id, content_id = action$content_id,
        expected_text = action$expected_text, plot_identity = action$content_id,
        input_id = ns(ui$ack), text_id = ns(ui$text),
        plot_id = if (is.null(ui$plot)) NULL else ns(ui$plot))
}

#' Bound the completion context: plain data only, no live owners, at most 4 MiB.
#' @noRd
validateWorkflowCompletionContext <- function(context) {
    if (is.null(context)) return(NULL)
    if (!is.list(context) || !is.list(context$peptide_qc) || hasMutableWorkflowEntityOwner(context) ||
        as.numeric(object.size(context)) > 4194304)
        cli::cli_abort("unsupported completion context", class = "workflow_execution_invalid")
    context
}

#' The workflow owner the controller finishes a worker result on: the project owner plus the GUI's
#' accumulated QC parameters and prior progress from the completion context (empty without one).
#' @noRd
newWorkflowCompletionWorkflow <- function(manager, context = NULL) {
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- manager
    workflow$qc_params <- list(peptide_qc = context$peptide_qc %||% list())
    if (is.list(context$prior_progress) && length(context$prior_progress))
        workflow$filtering_progress <- do.call(methods::new, c(list(Class = "FilteringProgress"), context$prior_progress))
    workflow
}

#' The progress product of a peptide QC state: `qc_progress` for the family and the sample filter,
#' `qvalue_progress` for the Q-value step (KIT-014).
#' @noRd
workflowPeptideProgressProduct <- function(record) {
    product <- record$products$qc_progress %||% record$products$qvalue_progress
    if (is.null(product)) cli::cli_abort("peptide progress product is missing", class = "workflow_execution_invalid")
    product
}

#' Recover the bounded display product of one accepted managed peptide action (the peptide operations
#' publish their progress product under `qc_progress`, the Q-value step under `qvalue_progress`).
#' @param manager Reference-backed owner.
#' @param result Accepted worker/controller result.
#' @param workflow_data Existing GUI workflow bus.
#' @return The display handoff (`stateReference`, text, summary, `qcParams`) without assay hydration.
#' @noRd
readManagedPeptideResult <- function(manager, result, workflow_data) {
    record <- readDIAStateMetadata(manager, result$state_reference)
    product <- workflowPeptideProgressProduct(record$data)
    envelope <- hydrateSavedWorkflowEntity(manager$.context$store, product)
    checkSavedWorkflowConsumerPayload("scientific_product", envelope)
    value <- c(list(stateReference = result$state_reference), envelope$value)
    # Merge the step's entries into the bus (audit PE-02): earlier steps' parameters are kept.
    qc_params <- shiny::isolate(workflow_data$qc_params) %||% list()
    if (is.null(qc_params$peptide_qc)) qc_params$peptide_qc <- list()
    for (name in names(value$qcParams$peptide_qc)) qc_params$peptide_qc[[name]] <- value$qcParams$peptide_qc[[name]]
    workflow_data$qc_params <- qc_params
    value
}

#' Install one module's managed-step state machine (KIT-008: the reference module pattern shared by the
#' precursor rollup and the sample filter). The step reserves the application's one heavy action, runs
#' the registered operation on the managed worker, polls it, delivers the accepted product to the module's
#' outputs, and closes the action only after the browser acknowledges the delivered content.
#' @param session,input,workflow_data The module's session, inputs and the workflow bus.
#' @param operation Registered operation name.
#' @param publication Immutable publication context (`experiment_label`, `experiment_paths`, `omic_type`).
#' @param ui Output names inside the module namespace: `text`, `plot`, the `ack` input, the working
#'   `notification` id, the `label` for messages and the `step_label` that identifies the content.
#' @param deliver Function(result) that renders the accepted result into the module's outputs and
#'   returns the plot grid.
#' @param read_result Function(manager, accepted, workflow_data) returning the display handoff.
#' @return An environment: `apply(parameters)` starts or retries an action; `active()` the current
#'   action or NULL; `cancel(reason)` releases a running action (Revert, session end).
#' @noRd
installWorkflowManagedStep <- function(session, input, workflow_data, operation, publication, ui, deliver,
    read_result = readManagedPeptideResult) {
    step <- new.env(parent = emptyenv())
    managed_action <- shiny::reactiveVal(NULL)
    pending_publication <- shiny::reactiveVal(NULL)
    working <- paste("Applying", ui$label, "...")
    failure_prefix <- paste0("Error applying ", ui$label, ":")
    ready_message <- function(action) session$sendCustomMessage("multischolar-perf107-ready",
        workflowManagedReadyPayload(session$ns, action, ui))
    fail <- function(action, error, status = "error") {
        if (identical(status, "publication_pending")) pending_publication(action)
        if (!is.null(action$owner) && !is.null(action$owner$job)) cancelWorkflowAction(action$owner, status)
        if (!is.null(action$reservation) && !action$reservation$closed) finishWorkflowGUIAction(action$reservation, status)
        managed_action(NULL)
        shiny::removeNotification(ui$notification)
        message <- paste(failure_prefix, conditionMessage(error))
        logger::log_error(message)
        shiny::showNotification(message, type = "error", duration = 15)
        invisible(NULL)
    }
    deliver_action <- function(action, result) {
        plot_grid <- withDIAStateScope(workflow_data$state_manager, function() deliver(result))
        content_id <- digest::digest(list(result$resultText,
            result$stateReference$entity_id %||% workflow_data$state_manager$current_state, ui$step_label), algo = "sha256")
        action$result <- result
        action$plot_grid <- plot_grid
        action$content_id <- content_id
        action$accepted_manager <- workflow_data$state_manager
        action$accepted_reference <- getDIAStateReference(workflow_data$state_manager)
        # Match renderText's output conversion exactly, including its treatment of a final newline.
        action$expected_text <- paste(utils::capture.output(cat(result$resultText)), collapse = "\n")
        if (is.function(ui$plot_identity)) ui$plot_identity(content_id)
        action$ack_started <- proc.time()[["elapsed"]]
        action$phase <- "awaiting_browser"
        managed_action(action)
        session$onFlushed(function() ready_message(action), once = TRUE)
        invisible(NULL)
    }
    step$active <- function() managed_action()
    step$pending <- function() pending_publication()
    step$cancel <- function(reason) {
        pending_publication(NULL)
        action <- shiny::isolate(managed_action())
        if (is.null(action)) return(invisible(NULL))
        if (!is.null(action$owner) && !is.null(action$owner$job)) cancelWorkflowAction(action$owner, reason)
        if (!action$reservation$closed) finishWorkflowGUIAction(action$reservation, reason)
        managed_action(NULL)
        shiny::removeNotification(ui$notification)
        invisible(NULL)
    }
    step$apply <- function(parameters) {
        application <- workflowSessionApplication(session)
        if (is.null(application)) cli::cli_abort("no managed application on this session", class = "workflow_execution_scope")
        if (!is.null(shiny::isolate(managed_action()))) {
            shiny::showNotification(paste("A", ui$label, "action is already running"), type = "warning")
            return(invisible(NULL))
        }
        action <- list(owner = NULL, phase = "starting")
        action$reservation <- tryCatch(beginWorkflowGUIAction(application, workflow_data$state_manager,
            operation, publication), error = function(error) error)
        if (inherits(action$reservation, "error")) {
            shiny::showNotification(conditionMessage(action$reservation), type = "error", duration = 15)
            return(invisible(NULL))
        }
        shiny::showNotification(working, id = ui$notification, duration = NULL)
        tryCatch({
            manager <- workflow_data$state_manager
            if (!inherits(manager, "DIAStateAdapter") || !usesDIAReferenceStates(manager))
                stop("the managed worker needs the on-disk state owner (DIA-NN import on Linux)")
            pending <- shiny::isolate(pending_publication())
            if (!is.null(pending)) {
                if (!identical(pending$owner$manager, manager) || !identical(pending$parameters, parameters))
                    stop("retry pending publication with its original parameters first")
                action$owner <- pending$owner
                action$worker_request_id <- pending$worker_request_id
                action$parameters <- pending$parameters
                accepted <- resumeWorkflowAction(action$owner, action$worker_request_id)
                deliver_action(action, read_result(manager, accepted, workflow_data))
                pending_publication(NULL)
                return(invisible(NULL))
            }
            action$parameters <- parameters
            owner <- manager$.context$store$execution_owner
            if (is.null(owner)) owner <- newWorkflowExecutionController(manager, application)
            action$owner <- owner
            action$phase <- "running"
            action$worker_request_id <- runWorkflowAction(operation, list(state = getDIAStateReference(manager)),
                parameters, owner, reservation = action$reservation,
                completion_context = workflowManagedCompletionContext(operation, workflow_data))
            managed_action(action)
        }, error = function(error) fail(action, error))
        invisible(NULL)
    }
    shiny::observe({
        action <- managed_action()
        if (is.null(action) || !identical(action$phase, "running")) return()
        if (!identical(workflow_data$state_manager, action$owner$manager)) {
            fail(action, simpleError("project was replaced"), "project_replaced")
            return()
        }
        shiny::invalidateLater(50, session)
        result <- tryCatch(pollWorkflowAction(action$owner), error = function(error) error)
        if (inherits(result, "error")) { fail(action, result); return() }
        if (result$status %in% c("admitting", "running")) return()
        if (!identical(result$status, "complete")) {
            fail(action, simpleError(paste("managed action ended:", result$status)), result$status)
            return()
        }
        tryCatch({
            accepted <- resumeWorkflowAction(action$owner, action$worker_request_id)
            deliver_action(action, read_result(workflow_data$state_manager, accepted, workflow_data))
        }, error = function(error) fail(action, error))
    })
    shiny::observe({
        action <- managed_action()
        if (is.null(action) || !identical(action$phase, "awaiting_browser")) return()
        if (!identical(workflow_data$state_manager, action$accepted_manager) ||
            !identical(getDIAStateReference(workflow_data$state_manager), action$accepted_reference)) {
            fail(action, simpleError("result was superseded"), "stale")
            return()
        }
        if (proc.time()[["elapsed"]] - action$ack_started >= 30) {
            fail(action, simpleError("browser acknowledgement timed out"), "browser_ack_timeout")
            return()
        }
        shiny::invalidateLater(1000, session)
        ready_message(action)
    })
    shiny::observeEvent(input[[ui$ack]], {
        action <- managed_action()
        acknowledgement <- input[[ui$ack]]
        if (is.null(action) || !identical(action$phase, "awaiting_browser") || !is.list(acknowledgement) ||
            !identical(acknowledgement$request_id, action$reservation$id) ||
            !identical(acknowledgement$content_id, action$content_id) || !isTRUE(acknowledgement$plot_ready) ||
            !is.numeric(acknowledgement$stable_ms) || length(acknowledgement$stable_ms) != 1L ||
            !is.finite(acknowledgement$stable_ms) || acknowledgement$stable_ms < 200) return()
        if (!identical(workflow_data$state_manager, action$accepted_manager) ||
            !identical(getDIAStateReference(workflow_data$state_manager), action$accepted_reference)) {
            fail(action, simpleError("result was superseded"), "stale")
            return()
        }
        receipt <- tryCatch(finishWorkflowGUIAction(action$reservation, "complete",
            list(request_id = acknowledgement$request_id, content_id = acknowledgement$content_id,
                plot_ready = TRUE, stable_seconds = acknowledgement$stable_ms / 1000)), error = function(error) error)
        if (inherits(receipt, "error")) { fail(action, receipt); return() }
        managed_action(NULL)
        shiny::removeNotification(ui$notification)
        logger::log_info(paste(ui$label, "applied successfully"))
        shiny::showNotification(paste(ui$label, "applied successfully"), type = "message")
    })
    if (is.function(session$onSessionEnded)) session$onSessionEnded(function() step$cancel("session_ended"))
    step
}

#' Open the one application-wide GUI action before reading scientific payloads.
#' @param application Managed application owner attached to the Shiny session.
#' @param manager Current native or reference-backed state owner.
#' @param operation Registered operation identifier.
#' @param publication Immutable real project/output context.
#' @return Process-local reservation borrowed by an optional worker action.
#' @noRd
beginWorkflowGUIAction <- function(application, manager, operation, publication) {
    entered <- proc.time()[["elapsed"]]
    checkWorkflowApplication(application)
    entry <- assertWorkflowOperationDispatchable(operation)
    if (!is.list(publication) || !identical(sort(names(publication)),
        sort(c("experiment_label", "experiment_paths", "omic_type"))) ||
        !is.character(publication$experiment_label) || length(publication$experiment_label) != 1L ||
        is.na(publication$experiment_label) || !nzchar(publication$experiment_label) ||
        !identical(publication$omic_type, entry$omic_type) ||
        !is.list(publication$experiment_paths) || hasMutableWorkflowEntityOwner(publication) ||
        as.numeric(object.size(publication)) > 65536) {
        cli::cli_abort("invalid GUI publication context", class = "workflow_execution_invalid")
    }
    if (!is.null(application$active)) {
        cli::cli_abort("another heavy action is active", class = "workflow_execution_busy")
    }
    peak <- openWorkflowPeakWindow(application)
    admission <- filelock::lock(file.path(application$config$root,
        ".heavy-action.lock"), timeout = 0)
    if (is.null(admission)) {
        close(peak)
        cli::cli_abort("application already has a heavy action", class = "workflow_execution_busy")
    }
    ready <- FALSE
    on.exit(if (!ready) {
        filelock::unlock(admission)
        close(peak)
    })
    accounting <- workflowCgroupAccounting(application$cgroup)
    id <- paste0("gui-", substr(digest::digest(list(Sys.time(), tempfile(),
        application$config$stem), algo = "sha256"), 1L, 32L))
    capture <- if (inherits(manager, "DIAStateAdapter")) {
        reference <- getDIAStateReference(manager)
        list(owner = manager$.context$store$identity, head = reference$entity_id,
            history = as.character(unlist(manager$getHistory())),
            current_state = manager$current_state)
    } else {
        list(owner = list(class = class(manager)), head = NULL,
            history = as.character(unlist(manager$getHistory())),
            current_state = manager$current_state)
    }
    reservation <- new.env(parent = emptyenv())
    reservation$schema <- 1L
    reservation$id <- id
    reservation$application <- application
    reservation$manager <- manager
    reservation$operation <- operation
    reservation$publication <- publication
    reservation$capture <- capture
    reservation$started <- entered
    reservation$peak <- peak
    reservation$accounting <- accounting
    reservation$admission <- admission
    reservation$pid <- Sys.getpid()
    reservation$closed <- FALSE
    reservation$worker_request_id <- NULL
    class(reservation) <- c("WorkflowGUIActionReservation", "environment")
    application$active <- id
    ready <- TRUE
    reservation
}

#' Close the GUI-complete action after current browser content is acknowledged.
#' @param reservation Process-local beginWorkflowGUIAction result.
#' @param status Terminal classification.
#' @param acknowledgement Bounded current request/head browser acknowledgement.
#' @return Durable, assay-free complete-operation receipt.
#' @noRd
finishWorkflowGUIAction <- function(reservation, status, acknowledgement = NULL) {
    if (!inherits(reservation, "WorkflowGUIActionReservation") || reservation$closed ||
        !identical(reservation$pid, Sys.getpid())) {
        cli::cli_abort("GUI action reservation is unavailable", class = "workflow_execution_stale")
    }
    application <- reservation$application
    if (!identical(application$active, reservation$id)) {
        cli::cli_abort("GUI action was superseded", class = "workflow_execution_stale")
    }
    success <- identical(status, "complete")
    valid_ack <- is.list(acknowledgement) && identical(acknowledgement$request_id,
        reservation$id) && is.character(acknowledgement$content_id) &&
        length(acknowledgement$content_id) == 1L && !is.na(acknowledgement$content_id) &&
        nzchar(acknowledgement$content_id) && isTRUE(acknowledgement$plot_ready) &&
        is.numeric(acknowledgement$stable_seconds) &&
        length(acknowledgement$stable_seconds) == 1L &&
        is.finite(acknowledgement$stable_seconds) && acknowledgement$stable_seconds >= 0.2
    if (success && !valid_ack) {
        cli::cli_abort("current browser acknowledgement required",
            class = "workflow_execution_incomplete")
    }
    peak_closed <- FALSE
    on.exit({
        if (!peak_closed) close(reservation$peak)
        if (!is.null(reservation$admission)) filelock::unlock(reservation$admission)
        reservation$admission <- NULL
        reservation$closed <- TRUE
        if (identical(application$active, reservation$id)) application$active <- NULL
    }, add = TRUE)
    browser <- workflowBrowserMembership(application)
    monitor <- NULL
    if (identical(application$config$mode, "gui")) {
        snapshot_id <- sub("^gui-", "wfx-", reservation$id)
        writeWorkflowExecutionMessage(file.path(application$config$root,
            "application-monitor-snapshot-request.rds"),
            list(schema = 1L, request_id = snapshot_id))
        snapshot_path <- file.path(application$config$root,
            paste0("application-monitor-snapshot-", snapshot_id, ".rds"))
        deadline <- proc.time()[["elapsed"]] + 2
        while (!file.exists(snapshot_path) && proc.time()[["elapsed"]] < deadline) Sys.sleep(0.02)
        if (file.exists(snapshot_path)) monitor <- readRDS(snapshot_path)
    }
    # Once finishing starts, IO/accounting failure must not strand the admission
    # lock or permit a second attempt to close the same peak descriptor.
    peak_closed <- TRUE
    peak <- closeWorkflowPeakWindow(reservation$peak)
    accounting <- workflowCgroupAccounting(application$cgroup)
    receipt <- list(schema = 1L, request_id = reservation$id, status = status,
        operation = reservation$operation,
        worker_request_id = reservation$worker_request_id,
        capture = reservation$capture,
        elapsed_seconds = proc.time()[["elapsed"]] - reservation$started,
        peak_bytes = peak, whole_service_peak_bytes = workflowMemoryCounter(file.path(
            application$cgroup, "memory.peak")), accounting_delta = Map(`-`, accounting,
            reservation$accounting), browser_acknowledged = valid_ack,
        browser_memory_included = isTRUE(browser$valid), browser = browser,
        monitor = monitor,
        monitor_basis = "application lifetime through final GUI acknowledgement; conservative attribution",
        measurement_valid = !is.null(monitor) && is.null(monitor$reason) &&
            is.numeric(monitor$max_observation_interval) &&
            length(monitor$max_observation_interval) == 1L &&
            is.finite(monitor$max_observation_interval) && monitor$max_observation_interval <= 0.1,
        acknowledgement = acknowledgement)
    path <- file.path(application$config$root, paste0(reservation$id, "-gui-completion.rds"))
    writeWorkflowExecutionMessage(path, receipt)
    filelock::unlock(reservation$admission)
    reservation$admission <- NULL
    reservation$closed <- TRUE
    application$active <- NULL
    receipt
}

#' Validate the proof driver's real service membership without inspecting payloads.
#' @param application Managed application owner.
#' @return Bounded membership result.
#' @noRd
workflowBrowserMembership <- function(application) {
    held <- list(valid = FALSE)
    path <- file.path(application$config$root, "application-browser-member.rds")
    if (!file.exists(path) || file.size(path) > 65536 || nzchar(Sys.readlink(path))) {
        return(held)
    }
    member <- readRDS(path)
    if (!is.list(member) || !identical(member$schema, 1L) ||
        !is.character(member$unit) || length(member$unit) != 1L ||
        !grepl("^multischolarwfx[a-f0-9]{16}-proof[.]service$", member$unit) ||
        !is.character(member$control_group) || length(member$control_group) != 1L ||
        !is.list(member$process)) return(held)
    actual <- workflowServiceProperties(member$unit)
    process <- tryCatch(workflowExecutionProcess(as.integer(member$process$pid)),
        error = function(error) NULL)
    valid <- isTRUE(actual$available) && actual$SubState %in% c("running", "start") &&
        identical(as.numeric(actual$MainPID), as.numeric(member$process$pid)) &&
        identical(actual$ControlGroup, member$control_group) &&
        startsWith(actual$ControlGroup, paste0(application$config$control_group, "/")) &&
        identical(process, member$process)
    list(valid = valid, unit = member$unit, control_group = member$control_group,
        process = member$process)
}

#' Check real cgroup membership and limits before project mutation or dispatch.
#' @param application Managed application resource owner.
#' @return TRUE invisibly.
#' @noRd
checkWorkflowApplication <- function(application) {
    if (!inherits(application, "WorkflowExecutionApplication") || application$closed ||
        !identical(application$pid, Sys.getpid()) || !identical(Sys.info()[["sysname"]], "Linux"))
        cli::cli_abort("managed application required", class = "workflow_execution_unsupported")
    config <- application$config
    actual <- workflowServiceProperties(config$controller_unit)
    monitor_ok <- TRUE
    if (identical(config$mode, "gui")) {
        monitor <- workflowServiceProperties(config$application_monitor_unit)
        monitor_ok <- isTRUE(monitor$available) && monitor$SubState %in% c("running", "start") &&
            startsWith(monitor$ControlGroup, paste0(config$control_group, "/"))
    }
    if (!isTRUE(actual$available) || !identical(as.numeric(actual$MainPID), as.numeric(Sys.getpid())) ||
        !startsWith(actual$ControlGroup, paste0(config$control_group, "/")) ||
        !monitor_ok ||
        workflowMemoryCounter(file.path(application$cgroup, "memory.max")) != config$policy$memory_max ||
        workflowMemoryCounter(file.path(application$cgroup, "memory.swap.max")) != 0)
        cli::cli_abort("application containment differs from policy", class = "workflow_execution_unsupported")
    invisible(TRUE)
}

#' Open a private action-window peak descriptor and reset it before invocation.
#' @param application Managed application owner.
#' @return Open descriptor; it is read once at action completion, then closed.
#' @noRd
openWorkflowPeakWindow <- function(application) {
    if (!inherits(application, "WorkflowExecutionApplication") || application$closed ||
        !identical(application$pid, Sys.getpid())) cli::cli_abort("managed peak owner required")
    refuse <- function(error) cli::cli_abort(paste("cannot reset the cgroup memory.peak counter;",
        "managed execution requires cgroup v2 on Linux 6.12 or later"),
        class = "workflow_execution_unsupported", parent = error)
    con <- tryCatch(suppressWarnings(file(file.path(application$cgroup, "memory.peak"), open = "r+")),
        error = refuse)
    ready <- FALSE
    on.exit(if (!ready) close(con))
    tryCatch(suppressWarnings({ writeLines("0", con); flush(con) }), error = refuse)
    ready <- TRUE
    con
}

#' Finish one FD-local peak window without subtracting GUI residency.
#' @param con Peak connection opened before action snapshot/admission.
#' @return Kernel-charged peak bytes.
#' @noRd
closeWorkflowPeakWindow <- function(con) {
    on.exit(close(con))
    seek(con, 0, origin = "start", rw = "read")
    value <- readLines(con, n = 1L, warn = FALSE)
    if (length(value) != 1L || !grepl("^[0-9]+$", value))
        cli::cli_abort("action peak counter unavailable", class = "workflow_execution_unsupported")
    as.numeric(value)
}

#' Start a managed Linux application with its aggregate cap set before its controller.
#' @param root Existing dedicated directory for bounded launch/journal metadata.
#' @param diagnostic Internal synthetic-lifecycle mode, never scientific qualification.
#' @param mode Either gui or one installed lifecycle suite.
#' @param test_root Source-bound test inputs for diagnostic mode only.
#' @param port Loopback Shiny port for the managed GUI.
#' @return Launch receipt with the exact owned slice/controller identities.
#' @noRd
launchWorkflowApplication <- function(root, diagnostic = FALSE, mode = "gui", test_root = NULL,
    port = 3838L, supervisor_unit = NULL) {
    if (!identical(Sys.info()[["sysname"]], "Linux") || !dir.exists(root) ||
        !mode %in% c("gui", "protocol", "resources", "adapters") || mode %in% c("protocol", "resources") && !isTRUE(diagnostic))
        cli::cli_abort("unsupported managed application launch", class = "workflow_execution_unsupported")
    policy <- workflowExecutionPolicy(diagnostic)
    if (workflowHostAvailable() < policy$startup_headroom)
        cli::cli_abort("insufficient host headroom", class = "workflow_execution_refused")
    workflowSystemd(c("show-environment"))
    if (!is.null(supervisor_unit)) {
        membership <- readLines("/proc/self/cgroup", warn = FALSE)
        if (!identical(mode, "gui") || !is.character(supervisor_unit) ||
            length(supervisor_unit) != 1L || is.na(supervisor_unit) ||
            !grepl("^multischolar-dia003-browser-smoke[.][A-Za-z0-9]{6}[.]service$", supervisor_unit) ||
            !any(endsWith(membership, paste0("/", supervisor_unit)))) {
            cli::cli_abort("supervisor is not this process's owned proof service",
                class = "workflow_execution_unsupported")
        }
    }
    id <- substr(digest::digest(list(Sys.time(), Sys.getpid(), tempfile()), algo = "sha256"), 1L, 16L)
    stem <- paste0("multischolarwfx", id)
    slice <- paste0(stem, ".slice")
    config <- list(schema = 1L, root = normalizePath(root), stem = stem, slice = slice,
        controller_unit = paste0(stem, "-controller.service"), policy = policy,
        build = workflowExecutionBuildIdentity(), library = dirname(system.file(package = "MultiScholaR")),
        mode = mode, port = as.integer(port), test_root = test_root)
    if (mode == "gui") {
        config$application_monitor_unit <- paste0(stem, "-appmonitor.service")
        member_path <- file.path(root, "gui-member.rds")
        if (file.exists(member_path)) {
            if (file.size(member_path) > 65536) {
                cli::cli_abort("GUI qualification member is oversized")
            }
            record <- workflowGuiLaunchRecord(readRDS(member_path),
                digest::digest(file = member_path, algo = "sha256"))
            config[names(record)] <- record
        }
    }
    if (mode != "gui") {
        test_file <- if (mode == "adapters") "test-workflow-execution-adapters.R" else "test-workflow-execution.R"
        config$test_file <- test_file
        inputs <- c(file.path(test_root, "tests/testthat", test_file),
            file.path(test_root, "tests/testthat/helper-workflow-checkpoint-fixtures.R"),
            file.path(test_root, "tests/testdata/workflow-execution/qualification.json"))
        if (any(!file.exists(inputs))) cli::cli_abort("installed lifecycle inputs are missing")
        if (mode == "adapters" && !diagnostic)
            inputs <- c(inputs, file.path(root, "timed-member.rds"),
                file.path(test_root, "tests/testthat/helper-workflow-bulk-rollup.R"))
        member_path <- file.path(root, "timed-member.rds")
        if (mode == "adapters" && file.exists(member_path)) {
            if (file.size(member_path) > 65536) cli::cli_abort("qualification member is oversized")
            member <- readRDS(member_path)
            if (identical(member$programme, "perf106-repair-v2")) {
                config$qualification <- list(programme = member$programme,
                    variant = member$variant,
                    member_sha256 = digest::digest(file = member_path, algo = "sha256"),
                    deadline = member$deadline)
                workflowQualificationVariant(config, member$operation)
                inputs <- unique(c(inputs, member_path))
            }
        }
        config$test_hashes <- stats::setNames(vapply(inputs, function(path)
            digest::digest(file = path, algo = "sha256"), character(1)), inputs)
    }
    if (identical(mode, "gui")) {
        # The outer supervisor must be able to clean up even if R is killed
        # between slice creation and publication of application.rds.
        writeLines(slice, file.path(root, "application-slice.txt"))
    }
    config$control_group <- createWorkflowApplicationSlice(stem, policy)
    ready <- FALSE
    on.exit(if (!ready) workflowSystemd(c("stop", slice), required = FALSE))
    group <- file.path("/sys/fs/cgroup", sub("^/", "", config$control_group))
    if (!nzchar(config$control_group) || workflowMemoryCounter(file.path(group, "memory.max")) != policy$memory_max ||
        workflowMemoryCounter(file.path(group, "memory.swap.max")) != 0)
        cli::cli_abort("aggregate memory controls are unavailable", class = "workflow_execution_unsupported")
    config_path <- file.path(root, "application.rds")
    if (file.exists(config_path)) cli::cli_abort("managed application directory already used")
    saveRDS(config, config_path, compress = FALSE, version = 3)
    controller_binds <- character()
    if (mode == "gui") {
        monitor_request <- list(schema = 1L, build = config$build, scope = "application",
            config = config, directory = config$root)
        monitor_path <- file.path(root, "application-monitor.rds")
        writeWorkflowExecutionMessage(monitor_path, monitor_request)
        workflowLaunchService(config, config$application_monitor_unit, "monitor", monitor_path,
            digest::digest(file = monitor_path, algo = "sha256"),
            binds = supervisor_unit %||% character(), deadline = 86500L)
        ready_path <- file.path(root, "application-monitor-ready.rds")
        deadline <- proc.time()[["elapsed"]] + 30
        while (!file.exists(ready_path) && proc.time()[["elapsed"]] < deadline) Sys.sleep(0.02)
        if (!file.exists(ready_path)) {
            cli::cli_abort("application monitor did not start", class = "workflow_execution_unsupported")
        }
        controller_binds <- config$application_monitor_unit
    }
    workflowLaunchService(config, config$controller_unit, "application", config_path,
        digest::digest(file = config_path, algo = "sha256"), deadline = if (!is.null(config$qualification))
            config$qualification$deadline else if (!is.null(config$journey)) config$journey$deadline else if (mode == "adapters") { if (diagnostic) 300L else 1800L } else
            if (mode == "resources" && diagnostic) 900L else if (diagnostic) 600L else 86400L,
        binds = controller_binds)
    ready <- TRUE
    config
}

#' Launch only package-defined service entry points beneath the owned application slice.
#' @param config Frozen managed application configuration.
#' @param unit Owned service name.
#' @param role Fixed internal entry role.
#' @param path,hash Bounded request file and expected hash.
#' @param cap Worker memory limit, or NULL for a controller/monitor.
#' @param binds Services whose loss must stop this service.
#' @param deadline Hard service lifetime in seconds.
#' @return NULL invisibly after launch acceptance, not a completion claim.
#' @noRd
workflowLaunchService <- function(config, unit, role, path, hash, cap = NULL, binds = character(), deadline = 300L) {
    if (!role %in% c("application", "monitor", "worker")) cli::cli_abort("invalid service role")
    expression <- switch(role,
        application = "MultiScholaR:::workflowManagedMain(commandArgs(TRUE))",
        monitor = "MultiScholaR:::workflowMonitorMain(commandArgs(TRUE))",
        worker = "MultiScholaR:::workflowWorkerMain(commandArgs(TRUE))")
    args <- c("--unit", unit, "--slice", config$slice,
        "--property=Type=exec", "--property=KillMode=control-group", "--property=TimeoutStopSec=1",
        "--property=MemoryAccounting=yes", "--property=CPUAccounting=yes", "--property=IOAccounting=yes",
        "--property=MemorySwapMax=0", paste0("--property=RuntimeMaxSec=", deadline),
        paste0("--property=StandardOutput=append:", file.path(config$root, paste0(role, "-", basename(path), ".log"))),
        "--property=StandardError=inherit", paste0("--setenv=R_LIBS=", paste(unique(c(config$library, .libPaths())), collapse = .Platform$path.sep)),
        "--setenv=OMP_NUM_THREADS=1", "--setenv=OPENBLAS_NUM_THREADS=1", "--setenv=MKL_NUM_THREADS=1",
        "--setenv=DUCKPLYR_FALLBACK_AUTOUPLOAD=0", "--setenv=DUCKPLYR_FALLBACK_COLLECT=0")
    if (role == "worker") args <- c(args, "--remain-after-exit")
    if (role == "application" && (config$mode %in% c("protocol", "resources") ||
        config$mode == "adapters" && isTRUE(config$policy$diagnostic))) args <- c(args, "--wait")
    if (!is.null(cap)) args <- c(args, paste0("--property=MemoryMax=", sprintf("%.0f", cap)), "--property=OOMPolicy=kill")
    if (length(binds)) args <- c(args, paste0("--property=BindsTo=", paste(binds, collapse = " ")),
        paste0("--property=After=", paste(binds, collapse = " ")))
    scratch <- file.path(dirname(path), paste0(role, "-scratch"))
    dir.create(scratch, showWarnings = FALSE)
    args <- c(args, paste0("--setenv=TMPDIR=", scratch))
    if (role == "worker") args <- c(args, "--setenv=MULTISCHOLAR_WORKFLOW_WORKER=1")
    args <- c(args, file.path(R.home("bin"), "Rscript"), "--vanilla", "-e", expression, path, hash)
    workflowSystemd(args, command = "systemd-run")
    invisible(NULL)
}

#' Read a hashed, bounded process request before interpreting it.
#' @param arguments Request filename/hash passed as separate argv values.
#' @return Validated metadata list.
#' @noRd
readWorkflowExecutionRequest <- function(arguments) {
    if (length(arguments) != 2L || !file.exists(arguments[[1L]]) || file.size(arguments[[1L]]) > 1048576 ||
        !identical(digest::digest(file = arguments[[1L]], algo = "sha256"), arguments[[2L]]))
        cli::cli_abort("execution request is missing or changed", class = "workflow_execution_invalid")
    request <- readRDS(arguments[[1L]])
    if (!is.list(request) || !identical(request$schema, 1L) ||
        !identical(request$build, workflowExecutionBuildIdentity()))
        cli::cli_abort("execution schema/build differs", class = "workflow_execution_invalid")
    request
}

#' Bootstrap the managed controller; diagnostic tests run against the installed namespace.
#' @param arguments Frozen application path/hash.
#' @return No value while the application is running.
#' @noRd
workflowManagedMain <- function(arguments) {
    config <- readWorkflowExecutionRequest(arguments)
    if (!identical(config$policy, workflowExecutionPolicy(config$policy$diagnostic)))
        cli::cli_abort("application policy differs from server policy")
    app <- WorkflowExecutionApplication$new()
    app$config <- config
    app$cgroup <- file.path("/sys/fs/cgroup", sub("^/", "", config$control_group))
    app$pid <- Sys.getpid()
    checkWorkflowApplication(app)
    on.exit(closeWorkflowManagedApplication(app))
    if (identical(config$mode, "fault_case")) {
        workflowControllerFaultMain(app)
        return(invisible(NULL))
    }
    if (config$mode != "gui") {
        for (path in names(config$test_hashes)) if (!identical(
            digest::digest(file = path, algo = "sha256"), config$test_hashes[[path]]))
            cli::cli_abort("lifecycle input changed")
        env <- new.env(parent = asNamespace("MultiScholaR"))
        env$workflow_test_application <- app
        env$workflow_test_mode <- config$mode
        env$workflow_test_root <- config$test_root
        sys.source(file.path(config$test_root, "tests/testthat/helper-workflow-checkpoint-fixtures.R"), env)
        results <- testthat::test_file(file.path(config$test_root, "tests/testthat",
            config$test_file %||% "test-workflow-execution.R"),
            env = env, load_helpers = FALSE, reporter = "summary", stop_on_failure = FALSE)
        results <- as.data.frame(results)
        saveRDS(results, file.path(config$root, "lifecycle-tests.rds"))
        if (any(results$failed > 0 | results$error)) stop("installed lifecycle checks failed")
        return(invisible(NULL))
    }
    # The native DIA session path is the production default (KIT-011); a KIT journey record no longer
    # has to switch it on (`config$journey$native_sessions` stays recorded for the harness).
    shiny::runApp(run_app(workflow_execution = app), host = "127.0.0.1",
        port = config$port, launch.browser = FALSE)
}

#' Interpret a bounded GUI launch record (`gui-member.rds`) into frozen configuration fields.
#' @param member Record list; `kit-journey-v1` keeps the journey's deadline (at most 4 hours, protocol
#'   v1.2) and guard-injection flag. The PERF-107 qualification record is retired (KIT-008): production
#'   dispatch no longer depends on a qualification, and its lifecycle scenarios were retired at KIT-013.
#' @param member_sha256 Hash of the record file.
#' @return Named list of configuration fields (`journey` plus the guard flag).
#' @noRd
workflowGuiLaunchRecord <- function(member, member_sha256) {
    if (!is.list(member) || !identical(member$schema, 1L)) {
        cli::cli_abort("GUI launch record is invalid")
    }
    if (!identical(member$programme, "kit-journey-v1")) {
        cli::cli_abort("GUI launch record names an unsupported programme", class = "workflow_execution_unsupported")
    }
    deadline <- member$deadline
    if (!is.numeric(deadline) || length(deadline) != 1L || !is.finite(deadline) ||
        deadline < 60 || deadline > 14400) {
        cli::cli_abort("KIT journey record is invalid")
    }
    list(journey = list(programme = "kit-journey-v1", native_sessions = TRUE,
        deadline = as.numeric(deadline), member_sha256 = member_sha256),
        allow_guard_injection = isTRUE(member$guard_stop_test))
}

#' Close every controller, then always release the application monitor.
#' @param app Managed application owner at process exit.
#' @return NULL invisibly; the first controller failure is re-signalled after release.
#' @noRd
closeWorkflowManagedApplication <- function(app) {
    failure <- NULL
    for (owner in app$controllers) {
        outcome <- tryCatch({ closeWorkflowExecutionController(owner); NULL }, error = \(error) error)
        failure <- failure %||% outcome
    }
    if (identical(app$config$mode, "gui")) {
        writeWorkflowExecutionMessage(file.path(app$config$root,
            "application-monitor-done.rds"), list(schema = 1L))
    }
    app$closed <- TRUE
    if (!is.null(failure)) stop(failure)
    invisible(NULL)
}

#' Read a guard classification; a malformed payload is itself the classification.
#' @param path Existing guard-stop file.
#' @return One stop reason written by a package monitor, or `"guard_payload_invalid"`.
#' @noRd
readWorkflowGuardReason <- function(path) {
    # The complete set written by workflowMonitorMain's request and application scopes.
    known <- c("memory_guard", "host_headroom", "nonworker_reservation", "invalid_worker_registration",
        "invalid_monitor_snapshot", "invalid_guard_injection", "injected_guard_stop")
    guard <- tryCatch(readRDS(path), error = \(error) NULL)
    reason <- if (is.list(guard)) guard$reason
    if (is.character(reason) && length(reason) == 1L && reason %in% known) reason
    else "guard_payload_invalid"
}

#' Persist bounded process handoffs without changing existing scientific stores.
#' @param path Owned request-directory file.
#' @param value Plain metadata/result reference list.
#' @return NULL invisibly.
#' @noRd
writeWorkflowExecutionMessage <- function(path, value) {
    if (hasMutableWorkflowEntityOwner(value) || as.numeric(object.size(value)) > 8388608)
        cli::cli_abort("execution handoff must be bounded immutable metadata")
    tmp <- tempfile("message-", tmpdir = dirname(path))
    on.exit(if (file.exists(tmp)) unlink(tmp))
    saveRDS(value, tmp, compress = FALSE, version = 3)
    if (!file.rename(tmp, path)) cli::cli_abort("execution message publication failed")
    invisible(NULL)
}

#' Monitor outside the worker; systemd binds worker lifetime to this process and the controller.
#' @param arguments Frozen monitor request path/hash.
#' @return NULL after worker completion or an explicit guard stop.
#' @noRd
workflowMonitorMain <- function(arguments) {
    request <- readWorkflowExecutionRequest(arguments)
    group <- file.path("/sys/fs/cgroup", sub("^/", "", request$config$control_group))
    if (identical(request$scope, "application")) {
        ready <- file.path(request$directory, "application-monitor-ready.rds")
        writeWorkflowExecutionMessage(ready, list(schema = 1L, pid = Sys.getpid(),
            process = workflowExecutionProcess(), control_group = request$config$control_group))
        max_anon <- max_file <- 0
        max_interval <- 0
        previous <- proc.time()[["elapsed"]]
        stop_reason <- NULL
        repeat {
            now <- proc.time()[["elapsed"]]
            max_interval <- max(max_interval, now - previous)
            previous <- now
            current <- workflowMemoryCounter(file.path(group, "memory.current"))
            host <- workflowHostAvailable()
            stat <- readLines(file.path(group, "memory.stat"), warn = FALSE)
            metric <- function(key) as.numeric(sub(paste0("^", key, " "), "",
                grep(paste0("^", key, " "), stat, value = TRUE)))
            max_anon <- max(max_anon, metric("anon"))
            max_file <- max(max_file, metric("file"))
            if (current >= request$config$policy$early_stop) stop_reason <- "memory_guard"
            if (host < request$config$policy$host_floor) stop_reason <- "host_headroom"
            active_path <- file.path(request$directory, "application-active-worker.rds")
            # The controller unlinks the registration when its worker completes; a read that lands in that
            # window sees no file (KIT-017 screen run: the monitor died and took the application with it).
            active <- if (file.exists(active_path)) tryCatch(readRDS(active_path), error = function(e) NULL) else NULL
            if (!is.null(active)) {
                valid_active <- is.list(active) && identical(active$schema, 1L) &&
                    is.character(active$request_id) && length(active$request_id) == 1L &&
                    grepl("^wfx-[a-f0-9]{32}$", active$request_id) &&
                    is.character(active$worker_unit) && length(active$worker_unit) == 1L &&
                    is.numeric(active$reservation) && length(active$reservation) == 1L &&
                    is.finite(active$reservation) && active$reservation >= 0
                if (!valid_active) {
                    stop_reason <- "invalid_worker_registration"
                } else {
                    worker_group <- file.path(group, active$worker_unit)
                    if (file.exists(file.path(worker_group, "memory.current"))) {
                        worker <- workflowMemoryCounter(file.path(worker_group, "memory.current"))
                        if (current - worker > active$reservation) {
                            stop_reason <- "nonworker_reservation"
                        }
                    }
                }
            }
            snapshot_path <- file.path(request$directory, "application-monitor-snapshot-request.rds")
            if (file.exists(snapshot_path)) {
                snapshot <- readRDS(snapshot_path)
                if (is.list(snapshot) && identical(snapshot$schema, 1L) &&
                    is.character(snapshot$request_id) && length(snapshot$request_id) == 1L &&
                    grepl("^wfx-[a-f0-9]{32}$", snapshot$request_id)) {
                    writeWorkflowExecutionMessage(file.path(request$directory, paste0(
                        "application-monitor-snapshot-", snapshot$request_id, ".rds")),
                    list(schema = 1L, scope = "application", max_anon = max_anon,
                        max_file = max_file, max_observation_interval = max_interval,
                        reason = stop_reason))
                } else {
                    stop_reason <- "invalid_monitor_snapshot"
                }
                unlink(snapshot_path)
            }
            injection <- file.path(request$directory, "application-guard-inject.rds")
            if (file.exists(injection)) {
                injected <- readRDS(injection)
                if (!(isTRUE(request$config$policy$diagnostic) ||
                    isTRUE(request$config$allow_guard_injection)) ||
                    !identical(injected, list(schema = 1L, reason = "injected_guard_stop"))) {
                    stop_reason <- "invalid_guard_injection"
                } else {
                    stop_reason <- "injected_guard_stop"
                }
            }
            done <- file.path(request$directory, "application-monitor-done.rds")
            if (!is.null(stop_reason) || file.exists(done)) break
            Sys.sleep(request$config$policy$monitor_interval)
        }
        result <- list(schema = 1L, scope = "application", max_anon = max_anon,
            max_file = max_file, max_observation_interval = max_interval,
            reason = stop_reason)
        writeWorkflowExecutionMessage(file.path(request$directory,
            "application-monitor-result.rds"), result)
        if (!is.null(stop_reason)) {
            writeWorkflowExecutionMessage(file.path(request$directory,
                "application-guard-stop.rds"), c(result, list(
                    current = current, host_available = host)))
            workflowSystemd(c("stop", request$config$slice), required = FALSE)
        }
        return(invisible(NULL))
    }
    if (!is.null(request$scope)) {
        cli::cli_abort("unknown monitor scope", class = "workflow_execution_invalid")
    }
    ready <- file.path(request$directory, "monitor-ready.rds")
    writeWorkflowExecutionMessage(ready, list(schema = 1L, pid = Sys.getpid()))
    max_anon <- max_file <- 0
    max_interval <- 0
    previous <- proc.time()[["elapsed"]]
    stop_reason <- NULL
    repeat {
        now <- proc.time()[["elapsed"]]
        max_interval <- max(max_interval, now - previous)
        previous <- now
        current <- workflowMemoryCounter(file.path(group, "memory.current"))
        host <- workflowHostAvailable()
        stat <- readLines(file.path(group, "memory.stat"), warn = FALSE)
        metric <- function(key) as.numeric(sub(paste0("^", key, " "), "", grep(paste0("^", key, " "), stat, value = TRUE)))
        max_anon <- max(max_anon, metric("anon"))
        max_file <- max(max_file, metric("file"))
        if (current >= request$config$policy$early_stop) stop_reason <- "memory_guard"
        if (host < request$config$policy$host_floor) stop_reason <- "host_headroom"
        worker_group <- file.path(group, request$worker_unit)
        admission <- file.path(request$directory, "admitted.rds")
        if (file.exists(admission)) {
            admitted <- readRDS(admission)
            if (file.exists(file.path(worker_group, "memory.current"))) {
                worker <- workflowMemoryCounter(file.path(worker_group, "memory.current"))
                if (current - worker > admitted$reservation) stop_reason <- "nonworker_reservation"
            }
        }
        if (!is.null(stop_reason)) {
            writeWorkflowExecutionMessage(file.path(request$directory, "guard-stop.rds"),
                list(schema = 1L, reason = stop_reason, current = current, host_available = host))
            workflowSystemd(c("stop", request$worker_unit), required = FALSE)
            break
        }
        if (file.exists(file.path(request$directory, "monitor-done.rds"))) break
        Sys.sleep(request$config$policy$monitor_interval)
    }
    writeWorkflowExecutionMessage(file.path(request$directory, "monitor-result.rds"), list(schema = 1L,
        max_anon = max_anon, max_file = max_file, max_observation_interval = max_interval,
        reason = stop_reason))
    invisible(NULL)
}

#' Execute a finite diagnostic adapter; production science adapters are added by WFX-004.
#' @param request Frozen request containing only input references and bounded parameters.
#' @return Existing reference/summary/checkpoint structures, never a hydrated assay.
#' @noRd
computeWorkflowFixture <- function(request) {
    if (!isTRUE(request$config$policy$diagnostic) || !isTRUE(workflowActionSpec(request$operation)$diagnostic_only))
        cli::cli_abort("operation has no qualified adapter", class = "workflow_execution_unsupported")
    store <- newWorkflowEntityStore(request$project_root, request$owner$project_id, request$owner$workflow_id)
    on.exit(closeWorkflowEntityStore(store))
    if (!identical(store$identity, request$owner)) cli::cli_abort("project incarnation changed")
    record <- hydrateSavedWorkflowCheckpoint(store, request$input_references$state)
    mode <- request$parameters$mode
    writeWorkflowExecutionMessage(file.path(request$directory, "worker-started.rds"),
        list(schema = 1L, mode = mode, pid = Sys.getpid()))
    if (mode == "forbidden_commit") {
        manager <- newDIAStateAdapter(store, reopen = TRUE)
        commitDIAReferenceState(manager, "forbidden", record$data, record$config, "worker must refuse")
    }
    if (mode == "fail") stop("intentional bounded worker failure")
    if (mode == "sleep") Sys.sleep(request$parameters$seconds)
    if (mode == "allocate") {
        # A real child/native allocation inside W, not an RSS estimate or mocked guard.
        child <- parallel::mcparallel({
            writeWorkflowExecutionMessage(file.path(request$directory, "allocation-started.rds"),
                list(schema = 1L, pid = Sys.getpid(), bytes = request$parameters$bytes,
                    cgroup = readLines("/proc/self/cgroup", warn = FALSE)))
            values <- rep.int(1, request$parameters$bytes / 8)
            Sys.sleep(3)
            sum(values)
        }, mc.set.seed = FALSE)
        result <- parallel::mccollect(child)
        if (inherits(result[[1L]], "try-error") || is.null(result[[1L]])) stop("allocation child failed")
    }
    list(checkpoint = record$data, config = record$config, state_name = workflowFixtureOperation,
        reference = record$data$components$peptide_data, dependencies = record$data$dependencies,
        summary = list(operation = "identity_fixture", input = request$input_references$state$entity_id))
}

#' Worker entry point; it cannot invoke controller head mutations.
#' @param arguments Frozen request path/hash.
#' @return NULL after bounded success/failure handoff; uncaught process loss is recorded by the controller.
#' @noRd
workflowWorkerMain <- function(arguments) {
    request <- readWorkflowExecutionRequest(arguments)
    if (!identical(Sys.getenv("MULTISCHOLAR_WORKFLOW_WORKER"), "1")) stop("worker role is missing")
    result <- tryCatch({
        handler <- workflowActionSpec(request$operation)$compute
        # The registry's compute role selects the worker's entry: the fixture, or the scientific action for
        # the peptide and protein roles (KIT-010), whose handler the entry names.
        value <- switch(handler, fixture = computeWorkflowFixture(request),
            peptide = , protein = computeWorkflowPeptideAction(request),
            cli::cli_abort("operation adapter unavailable", class = "workflow_execution_unsupported"))
        list(schema = 1L, status = "complete", token = request$token, result = value)
    }, error = function(error) list(schema = 1L, status = "failed", token = request$token,
        error = substr(conditionMessage(error), 1L, 1024L)))
    writeWorkflowExecutionMessage(file.path(request$directory, "worker-result.rds"), result)
    if (result$status != "complete") stop(result$error, call. = FALSE)
    invisible(NULL)
}

#' Submit one finite heavy action using frozen references and the server resource owner.
#' @param operation Registered operation ID, never a function/executable.
#' @param input_references Named immutable descriptors.
#' @param parameters Small immutable operation parameters.
#' @param resource_policy Project controller created inside a managed application.
#' @param reservation Optional matching GUI action opened at handler entry.
#' @return Request ID immediately after monitor launch; pollWorkflowAction advances it.
#' @noRd
runWorkflowAction <- function(operation, input_references, parameters, resource_policy,
    reservation = NULL, completion_context = NULL) {
    entered <- proc.time()[["elapsed"]]
    owner <- resource_policy
    if (!inherits(owner, "WorkflowExecutionController") || owner$closed) cli::cli_abort("controller required")
    app <- owner$application
    outer <- !is.null(reservation)
    if (outer && (!inherits(reservation, "WorkflowGUIActionReservation") ||
        reservation$closed || !identical(reservation$pid, Sys.getpid()) ||
        !identical(reservation$application, app) || !identical(reservation$manager, owner$manager) ||
        !identical(reservation$operation, operation) ||
        !identical(app$active, reservation$id))) {
        cli::cli_abort("GUI action reservation differs", class = "workflow_execution_stale")
    }
    peak <- if (outer) reservation$peak else openWorkflowPeakWindow(app)
    transferred <- FALSE
    on.exit(if (!transferred && !outer) close(peak), add = TRUE)
    checkWorkflowApplication(app)
    assertDIAExportReady(owner$manager)
    if ((!outer && !is.null(app$active)) || !is.null(owner$job)) {
        cli::cli_abort("another heavy action is active", class = "workflow_execution_busy")
    }
    assertWorkflowOperationDispatchable(operation)
    validateWorkflowActionParameters(operation, parameters, app$config$policy, resolved = FALSE)
    if (!is.list(input_references) || !identical(names(input_references), "state") ||
        !isSavedWorkflowReference(input_references$state))
        cli::cli_abort("unsupported input reference", class = "workflow_execution_invalid")
    if (workflowHostAvailable() < app$config$policy$startup_headroom)
        cli::cli_abort("insufficient host headroom", class = "workflow_execution_refused")
    admission <- if (outer) reservation$admission else filelock::lock(file.path(
        app$config$root, ".heavy-action.lock"), timeout = 0)
    if (is.null(admission)) cli::cli_abort("application already has a heavy action",
        class = "workflow_execution_busy")
    ready <- FALSE
    on.exit(if (!ready) {
        if (!outer) filelock::unlock(admission)
        owner$job <- NULL
        if (!outer) app$active <- NULL
    }, add = TRUE)
    job <- new.env(parent = emptyenv())
    job$started <- entered
    job$peak <- peak
    job$outer_reservation <- if (outer) reservation else NULL
    job$accounting <- workflowCgroupAccounting(app$cgroup)
    store <- owner$manager$.context$store
    job$existing_entities <- list.files(file.path(store$directory, "entities"))
    job$existing_pending <- list.files(file.path(store$directory, ".pending"), all.files = TRUE, no.. = TRUE)
    validateSavedWorkflowReference(store, input_references$state, payload = FALSE)
    if (!identical(input_references$state, getDIAStateReference(owner$manager)))
        cli::cli_abort("action input is not the selected immutable checkpoint", class = "workflow_execution_stale")
    validateWorkflowActionInput(operation, owner$manager, input_references)
    parameters <- resolveWorkflowActionParameters(operation, owner$manager,
        input_references$state, parameters, app$config$policy)
    id <- paste0("wfx-", substr(digest::digest(list(Sys.time(), tempfile(), owner$control$generation), algo = "sha256"), 1, 32))
    paths <- workflowExportPaths(store$directory, id)
    token <- list(request_id = id, incarnation = store$identity$incarnation,
        input_ids = vapply(input_references, `[[`, character(1), "entity_id"),
        expected_head = workflowExecutionHead(store)$entity_id, generation = owner$control$generation,
        controller_id = owner$control$controller_id, parameters = digest::digest(parameters, algo = "sha256"),
        execution = workflowExecutionBuildIdentity())
    variant <- workflowQualificationVariant(app$config, operation)
    if (!is.null(app$config$qualification)) token$qualification <- digest::digest(
        app$config$qualification, algo = "sha256")
    request <- list(schema = 1L, build = token$execution, operation = operation, input_references = input_references,
        parameters = parameters, token = token, owner = store$identity, project_root = dirname(store$directory),
        directory = paths$request_dir, config = app$config)
    job$request <- request
    # KIT-008 (audit PE-02): the GUI's accumulated QC parameters and prior progress travel with the job so the
    # controller's completion writes the same progress product the in-process path writes.
    job$completion_context <- validateWorkflowCompletionContext(completion_context)
    job$paths <- paths
    job$admission <- admission
    job$stage <- "admitting"
    job$worker_unit <- paste0(app$config$stem, "-w", substr(token$request_id, 5, 16), ".service")
    job$monitor_borrowed <- identical(app$config$mode, "gui")
    job$monitor_unit <- if (job$monitor_borrowed) app$config$application_monitor_unit else
        paste0(app$config$stem, "-m", substr(token$request_id, 5, 16), ".service")
    job$worker_request <- file.path(paths$request_dir, "request.rds")
    writeWorkflowExecutionMessage(job$worker_request, request)
    if (outer) {
        writeWorkflowExecutionMessage(file.path(paths$request_dir,
            "publication-context.rds"), reservation$publication)
        reservation$worker_request_id <- id
    }
    monitor_path <- NULL
    if (!job$monitor_borrowed) {
        monitor <- list(schema = 1L, build = request$build, config = app$config,
            directory = paths$request_dir, worker_unit = job$worker_unit)
        monitor_path <- file.path(paths$request_dir, "monitor.rds")
        writeWorkflowExecutionMessage(monitor_path, monitor)
    }
    writeWorkflowExportJournal(paths, list(schema = "multischolar.workflow_action.v1", request_id = id,
        status = "admitting", expected_head = token$expected_head, token_hash = digest::digest(token, algo = "sha256")))
    owner$control$active_request <- id
    writeWorkflowExecutionControl(store, owner$control)
    owner$job <- job
    if (!outer) app$active <- id
    # Transfer ownership before launch: a failed acknowledgement may still have created M.
    ready <- TRUE
    transferred <- TRUE
    launched <- tryCatch({
        if (!job$monitor_borrowed) {
            workflowLaunchService(app$config, job$monitor_unit, "monitor", monitor_path,
                digest::digest(file = monitor_path, algo = "sha256"), binds = app$config$controller_unit,
                deadline = app$config$policy$deadline + 60L)
        }
        TRUE
    }, workflow_execution_unsupported = function(error) {
        job$error <- substr(conditionMessage(error), 1L, 1024L)
        FALSE
    })
    if (!launched) {
        result <- finishWorkflowExecutionJob(owner, "admission_refused")
        cli::cli_abort("monitor launch failed: {job$error}",
            class = if (result$status == "cleanup_failed") "workflow_execution_cleanup_failed"
                else "workflow_execution_refused", request_id = id)
    }
    id
}

#' Poll a managed action without blocking the Shiny event loop on computation.
#' @param owner Project controller.
#' @return Bounded stage or completion/failure receipt.
#' @noRd
pollWorkflowAction <- function(owner) {
    if (is.null(owner$job)) return(list(status = "idle"))
    job <- owner$job
    app <- owner$application
    if (identical(job$stage, "cleanup_failed"))
        return(finishWorkflowExecutionJob(owner, job$terminal_status))
    if (!workflowExecutionStoreValid(owner$manager$.context$store))
        return(finishWorkflowExecutionJob(owner, "project_replaced"))
    guards <- c(file.path(job$paths$request_dir, "guard-stop.rds"),
        if (job$monitor_borrowed) file.path(app$config$root, "application-guard-stop.rds"))
    guard <- guards[file.exists(guards)][1L]
    if (length(guard) && !is.na(guard)) {
        return(finishWorkflowExecutionJob(owner, readWorkflowGuardReason(guard)))
    }
    result_path <- file.path(job$paths$request_dir, "worker-result.rds")
    monitor <- workflowServiceProperties(job$monitor_unit)
    if (!isTRUE(monitor$available) || !monitor$SubState %in% c("running", "start")) {
        # The worker's result is kept; only the measurement is lost with the monitor.
        if (!file.exists(result_path)) return(finishWorkflowExecutionJob(owner, "monitor_lost"))
        job$monitor_lost_after_result <- TRUE
    }
    if (job$stage == "admitting") {
        monitor_ready <- if (job$monitor_borrowed) file.path(app$config$root,
            "application-monitor-ready.rds") else file.path(job$paths$request_dir, "monitor-ready.rds")
        if (!file.exists(monitor_ready)) {
            if (proc.time()[["elapsed"]] - job$started > 30) return(finishWorkflowExecutionJob(owner, "monitor_start_timeout"))
            return(list(status = "admitting"))
        }
        admitted <- workflowAdmissionCap(app$config$policy,
            workflowMemoryCounter(file.path(app$cgroup, "memory.current")), workflowHostAvailable())
        reservation <- admitted$reservation
        cap <- admitted$cap
        # This finite fixture estimate includes R/package startup and bounded reader overhead.
        if (cap < 768 * 1024^2) return(finishWorkflowExecutionJob(owner, "admission_refused"))
        writeWorkflowExecutionMessage(file.path(job$paths$request_dir, "admitted.rds"),
            list(schema = 1L, reservation = reservation, worker_cap = cap, worker_cap_source = admitted$source))
        if (job$monitor_borrowed) {
            writeWorkflowExecutionMessage(file.path(app$config$root,
                "application-active-worker.rds"), list(schema = 1L,
                    request_id = job$request$token$request_id,
                    worker_unit = job$worker_unit, reservation = reservation))
        }
        job$worker_cap <- cap
        job$request$worker_cap <- cap
        writeWorkflowExecutionMessage(job$worker_request, job$request)
        job$stage <- "running"
        launched <- tryCatch({
            workflowLaunchService(app$config, job$worker_unit, "worker", job$worker_request,
                digest::digest(file = job$worker_request, algo = "sha256"), cap = cap,
                binds = c(app$config$controller_unit, job$monitor_unit), deadline = app$config$policy$deadline)
            TRUE
        }, workflow_execution_unsupported = function(error) {
            job$error <- substr(conditionMessage(error), 1L, 1024L)
            FALSE
        })
        if (!launched) return(finishWorkflowExecutionJob(owner, "admission_refused"))
        return(list(status = "running"))
    }
    state <- workflowServiceProperties(job$worker_unit)
    if (isTRUE(state$available) && state$SubState %in% c("running", "start")) return(list(status = "running"))
    if (!file.exists(result_path) && !any(file.exists(guards))) {
        # A worker may become terminal between two 50 ms monitor samples.
        # Give that already-running monitor two intervals to publish its
        # durable classification before treating the stop as an unguarded
        # worker failure.
        Sys.sleep(min(0.25, max(0.01, 2 * app$config$policy$monitor_interval)))
        guard <- guards[file.exists(guards)][1L]
        if (length(guard) && !is.na(guard)) {
            return(finishWorkflowExecutionJob(owner, readWorkflowGuardReason(guard)))
        }
        monitor <- workflowServiceProperties(job$monitor_unit)
        if (!isTRUE(monitor$available) || !monitor$SubState %in% c("running", "start"))
            return(finishWorkflowExecutionJob(owner, "monitor_lost"))
    }
    if (!isTRUE(state$available) || !identical(state$Result, "success") || !file.exists(result_path))
        return(finishWorkflowExecutionJob(owner, if (identical(state$Result, "oom-kill")) "memory_guard" else "worker_failed"))
    if (file.size(result_path) > 8388608) return(finishWorkflowExecutionJob(owner, "invalid_result"))
    result <- readRDS(result_path)
    if (!is.list(result) || !identical(result$schema, 1L) || !identical(result$status, "complete"))
        return(finishWorkflowExecutionJob(owner, "worker_failed"))
    acceptWorkflowActionResult(owner, result)
}

#' Accept exactly the frozen token inside the controller lease, then publish separately.
#' @param owner Leased controller.
#' @param message Immutable worker handoff.
#' @return Bounded terminal receipt.
#' @noRd
acceptWorkflowActionResult <- function(owner, message) {
    if (!inherits(owner, "WorkflowExecutionController") || owner$closed ||
        !identical(owner$pid, Sys.getpid()) || is.null(owner$lease))
        cli::cli_abort("completion controller is no longer active", class = "workflow_execution_stale")
    if (!workflowExecutionStoreValid(owner$manager$.context$store)) {
        if (!is.null(owner$job)) return(finishWorkflowExecutionJob(owner, "project_replaced"))
        cli::cli_abort("completion project was replaced", class = "workflow_execution_stale")
    }
    job <- owner$job
    if (!is.list(message) || !is.list(message$token) ||
        hasMutableWorkflowEntityOwner(message$token) || as.numeric(object.size(message$token)) > 65536)
        cli::cli_abort("invalid completion token", class = "workflow_execution_invalid")
    if (is.null(job) || !identical(message$token$request_id, job$request$token$request_id)) {
        id <- message$token$request_id
        if (!is.character(id) || length(id) != 1L || !grepl("^wfx-[a-f0-9]{32}$", id))
            cli::cli_abort("unknown completion", class = "workflow_execution_stale")
        if (!dir.exists(file.path(owner$manager$.context$store$directory, ".workflow-export", id)))
            cli::cli_abort("unknown completion", class = "workflow_execution_stale")
        paths <- workflowExportPaths(owner$manager$.context$store$directory, id)
        journal <- readWorkflowExecutionJournal(paths)
        if (is.null(journal) || !identical(journal$token_hash, digest::digest(message$token, algo = "sha256")) ||
            !journal$status %in% c("science_committed", "publication_pending", "complete"))
            cli::cli_abort("unknown completion", class = "workflow_execution_stale")
        path <- file.path(paths$request_dir, "accepted-result.rds")
        if (!file.exists(path) || !identical(digest::digest(file = path, algo = "sha256"), journal$result_sha256))
            cli::cli_abort("accepted completion changed", class = "workflow_execution_invalid")
        result <- readRDS(path)
        return(list(status = if (journal$status == "complete") "duplicate" else "publication_pending",
            result = result, delivery_current = is.null(job) &&
                identical(owner$control$generation, journal$delivery_generation) &&
                identical(getDIAStateReference(owner$manager), result$state_reference)))
    }
    store <- owner$manager$.context$store
    token <- job$request$token
    current <- workflowExecutionHead(store)
    if (!identical(message$token, token) || !identical(token$generation, owner$control$generation) ||
        !identical(token$controller_id, owner$control$controller_id) ||
        !identical(token$incarnation, readSavedStoreIdentity(file.path(store$directory, "store.json"))$incarnation) ||
        !identical(current$entity_id, token$expected_head)) return(finishWorkflowExecutionJob(owner, "stale"))
    value <- message$result
    if (hasMutableWorkflowEntityOwner(value) || !isDIAComponentCheckpoint(value$checkpoint) ||
        !isSavedWorkflowReference(value$reference) || as.numeric(object.size(value)) > 8388608)
        return(finishWorkflowExecutionJob(owner, "invalid_result"))
    validateSavedWorkflowReference(store, value$reference, payload = TRUE)
    writeWorkflowExecutionMessage(file.path(job$paths$request_dir, "pending-result.rds"), value)
    owner$accepting <- token$request_id
    on.exit(owner$accepting <- NULL)
    completeWorkflowAction(owner, job$request$operation, value, context = job$completion_context,
        variant = workflowQualificationVariant(job$request$config, job$request$operation), values = job$request$parameters)
    owner$accepting <- NULL
    journal <- readWorkflowExecutionJournal(job$paths)
    journal$status <- "publication_pending"
    writeWorkflowExportJournal(job$paths, journal)
    published <- tryCatch({ publishWorkflowFixtureResult(owner, job$paths); TRUE },
        error = function(error) {
            job$error <- substr(conditionMessage(error), 1L, 1024L)
            FALSE
        })
    finishWorkflowExecutionJob(owner, if (published) "complete" else "publication_pending")
}

#' Publish the fixture's bounded required output using the established output journal.
#' @param owner Leased controller.
#' @param paths Existing execution request paths.
#' @return Saved reference/summary handoff.
#' @noRd
publishWorkflowFixtureResult <- function(owner, paths) {
    journal <- readWorkflowExecutionJournal(paths)
    result_path <- file.path(paths$request_dir, "accepted-result.rds")
    if (!file.exists(result_path) || !identical(digest::digest(file = result_path, algo = "sha256"), journal$result_sha256))
        cli::cli_abort("accepted result receipt differs")
    result <- readRDS(result_path)
    if (!identical(getDIAStateReference(owner$manager), result$state_reference))
        cli::cli_abort("publication belongs to a superseded state", class = "workflow_execution_stale")
    publication <- workflowActionSpec(result$operation %||% workflowFixtureOperation)$publication
    if (publication %in% c("peptide", "protein", "da")) {
        if (publication == "da") publishWorkflowDaResult(owner, paths, result) else
        if (publication == "protein") publishWorkflowProteinResult(owner, paths, result) else
            publishWorkflowPeptideResult(owner, paths, result)
        journal$status <- "complete"
        journal$delivery_generation <- owner$control$generation
        writeWorkflowExportJournal(paths, journal)
        return(result)
    }
    output <- file.path(paths$request_dir, "required-output")
    dir.create(output, showWarnings = FALSE)
    owner$manager$.context$experiment_paths <- list(source_dir = output)
    publishDIAActionOutputs(owner$manager, result$state_reference, "wfx-fixture", list(output = output),
        function(directories) {
            workflowExecutionBoundary("publication")
            jsonlite::write_json(result$summary, file.path(directories$output, "summary.json"), auto_unbox = TRUE)
            list(status = "published")
        })
    journal$status <- "complete"
    journal$delivery_generation <- owner$control$generation
    writeWorkflowExportJournal(paths, journal)
    result
}

#' Retry only accepted publication; never dispatch or recompute scientific work.
#' @param owner Leased controller.
#' @param request_id Previously accepted idempotency token.
#' @return Original bounded result.
#' @noRd
resumeWorkflowAction <- function(owner, request_id) {
    if (!inherits(owner, "WorkflowExecutionController") || owner$closed ||
        !identical(owner$pid, Sys.getpid()) || is.null(owner$lease) ||
        !workflowExecutionStoreValid(owner$manager$.context$store))
        cli::cli_abort("publication controller is no longer active", class = "workflow_execution_stale")
    paths <- workflowExportPaths(owner$manager$.context$store$directory, request_id)
    journal <- readWorkflowExecutionJournal(paths)
    if (!journal$status %in% c("science_committed", "publication_pending", "complete"))
        cli::cli_abort("request has no accepted scientific result", class = "workflow_execution_invalid")
    if (journal$status != "complete") return(publishWorkflowFixtureResult(owner, paths))
    path <- file.path(paths$request_dir, "accepted-result.rds")
    if (!identical(digest::digest(file = path, algo = "sha256"), journal$result_sha256)) cli::cli_abort("result changed")
    result <- readRDS(path)
    if (!identical(getDIAStateReference(owner$manager), result$state_reference))
        cli::cli_abort("publication state was superseded", class = "workflow_execution_stale")
    journal$delivery_generation <- owner$control$generation
    writeWorkflowExportJournal(paths, journal)
    result
}

#' Fixed failure-injection boundary used only by owned lifecycle fixtures.
#' @param stage Bounded phase name.
#' @return NULL invisibly.
#' @noRd
workflowExecutionBoundary <- function(stage) {
    hook <- getOption("multischolar.workflow_execution_boundary")
    if (is.function(hook)) hook(stage)
    invisible(NULL)
}

#' Stop owned worker/monitor services and close one complete measurement window.
#' @param owner Leased controller.
#' @param status Explicit terminal classification; guards are failures, never success.
#' @return Bounded reference-free execution receipt.
#' @noRd
finishWorkflowExecutionJob <- function(owner, status) {
    job <- owner$job
    if (is.null(job)) return(list(status = "idle"))
    valid_store <- workflowExecutionStoreValid(owner$manager$.context$store)
    # Stop M first when its request paths no longer belong to this project.
    if (!valid_store && !job$monitor_borrowed) {
        workflowSystemd(c("stop", job$monitor_unit), required = FALSE)
    }
    worker <- workflowServiceProperties(job$worker_unit)
    workflowSystemd(c("stop", job$worker_unit), required = FALSE)
    if (valid_store && dir.exists(job$paths$request_dir) && !job$monitor_borrowed) {
        writeWorkflowExecutionMessage(file.path(job$paths$request_dir, "monitor-done.rds"), list(schema = 1L))
    }
    if (valid_store && job$monitor_borrowed) {
        writeWorkflowExecutionMessage(file.path(owner$application$config$root,
            "application-monitor-snapshot-request.rds"), list(schema = 1L,
                request_id = job$request$token$request_id))
    }
    deadline <- proc.time()[["elapsed"]] + 2
    monitor_path <- if (job$monitor_borrowed) file.path(owner$application$config$root,
        paste0("application-monitor-snapshot-", job$request$token$request_id, ".rds")) else
        file.path(job$paths$request_dir, "monitor-result.rds")
    while (valid_store && !file.exists(monitor_path) && proc.time()[["elapsed"]] < deadline) Sys.sleep(0.02)
    if (!job$monitor_borrowed) {
        workflowSystemd(c("stop", job$monitor_unit), required = FALSE)
    }
    worker_empty <- workflowServiceEmpty(owner$application, job$worker_unit)
    monitor_empty <- if (job$monitor_borrowed) TRUE else
        workflowServiceEmpty(owner$application, job$monitor_unit)
    if (!worker_empty || !monitor_empty) {
        job$stage <- "cleanup_failed"
        job$terminal_status <- status
        return(list(status = "cleanup_failed", request_id = job$request$token$request_id))
    }
    on.exit({
        if (is.null(job$outer_reservation)) filelock::unlock(job$admission)
        if (job$monitor_borrowed) {
            active_path <- file.path(owner$application$config$root,
                "application-active-worker.rds")
            if (file.exists(active_path)) {
                active <- readRDS(active_path)
                if (identical(active$request_id, job$request$token$request_id)) unlink(active_path)
            }
        }
        if (is.null(job$outer_reservation)) owner$application$active <- NULL
        owner$job <- NULL
    })
    roles <- if (job$monitor_borrowed) "worker" else c("worker", "monitor")
    if (valid_store) for (role in roles) {
        scratch <- file.path(job$paths$request_dir, paste0(role, "-scratch"))
        if (dir.exists(scratch) && !nzchar(Sys.readlink(scratch))) unlink(scratch, recursive = TRUE)
    }
    new_entities <- if (valid_store) setdiff(list.files(file.path(owner$manager$.context$store$directory,
        "entities")), job$existing_entities) else character()
    new_pending <- if (valid_store) setdiff(list.files(file.path(owner$manager$.context$store$directory,
        ".pending"), all.files = TRUE, no.. = TRUE), job$existing_pending) else character()
    artifact_roots <- if (valid_store) c(job$paths$request_dir,
        file.path(owner$manager$.context$store$directory, "entities", new_entities),
        file.path(owner$manager$.context$store$directory, ".pending", new_pending)) else character()
    files <- unique(unlist(lapply(artifact_roots, list.files,
        recursive = TRUE, full.names = TRUE, all.files = TRUE)))
    staged_bytes <- if (valid_store) sum(file.info(files)$size) else NA_real_
    peak <- if (is.null(job$outer_reservation)) closeWorkflowPeakWindow(job$peak) else NA_real_
    accounting <- workflowCgroupAccounting(owner$application$cgroup)
    delta <- Map(`-`, accounting, job$accounting)
    monitor <- if (valid_store && file.exists(monitor_path)) readRDS(monitor_path) else NULL
    receipt <- list(schema = 1L, request_id = job$request$token$request_id, status = status,
        token = job$request$token, error = job$error, elapsed_seconds = proc.time()[["elapsed"]] - job$started,
        peak_bytes = peak, whole_service_peak_bytes = workflowMemoryCounter(file.path(owner$application$cgroup, "memory.peak")),
        monitor = monitor, worker = worker, threads = 8L, browser_memory_included = FALSE,
        monitor_borrowed = job$monitor_borrowed,
        accounting_delta = delta, staged_bytes_at_teardown = staged_bytes,
        staged_basis = "new entity/map/summary payloads plus request outputs and owned pending files; retained input links excluded",
        new_entity_count = length(new_entities),
        spill_attempts = if (isTRUE(workflowActionSpec(job$request$operation)$diagnostic_only)) 0L else NA_integer_,
        spill_basis = if (isTRUE(workflowActionSpec(job$request$operation)$diagnostic_only))
            "identity fixture uses no database or spill-capable algorithm" else
            "engine spill disabled; attempted spill is not independently counted",
        measurement_valid = !is.null(monitor) && monitor$max_observation_interval <= 0.1 &&
            !isTRUE(job$monitor_lost_after_result),
        complete_operation_measurement = is.null(job$outer_reservation),
        cleanup = list(worker_stopped = worker_empty,
            monitor_stopped = if (job$monitor_borrowed) NA else monitor_empty,
            application_monitor_active = if (job$monitor_borrowed)
                isTRUE(workflowServiceProperties(job$monitor_unit)$available) else NA))
    if (isTRUE(job$monitor_lost_after_result)) receipt$monitor_lost_after_result <- TRUE
    valid_store <- workflowExecutionStoreValid(owner$manager$.context$store)
    if (!valid_store) {
        receipt$status <- "project_replaced"
        writeWorkflowExecutionMessage(file.path(owner$application$config$root,
            paste0(job$request$token$request_id, "-completion.rds")), receipt)
        return(receipt)
    }
    writeWorkflowExecutionMessage(file.path(job$paths$request_dir, "completion.rds"), receipt)
    journal <- readWorkflowExecutionJournal(job$paths)
    if (!journal$status %in% c("science_committed", "publication_pending", "complete")) journal$status <- status
    journal$completion_status <- status
    writeWorkflowExportJournal(job$paths, journal)
    owner$control$active_request <- NULL
    writeWorkflowExecutionControl(owner$manager$.context$store, owner$control)
    receipt
}

#' Bind only the managed application owner into the ordinary Shiny server.
#' @param application Managed application resource owner.
#' @return Existing server with session lifecycle ownership attached.
#' @noRd
workflowManagedServer <- function(application) {
    force(application)
    function(input, output, session) {
        session$userData$workflow_execution <- application
        app_server(input, output, session)
    }
}

#' Read aggregate CPU/I/O counters without overlapping process RSS arithmetic.
#' @param group Application cgroup path.
#' @return Kernel cumulative counters for a later delta.
#' @noRd
workflowCgroupAccounting <- function(group) {
    cpu <- readLines(file.path(group, "cpu.stat"), warn = FALSE)
    usage <- as.numeric(sub("^usage_usec ", "", grep("^usage_usec ", cpu, value = TRUE)))
    io <- readLines(file.path(group, "io.stat"), warn = FALSE)
    bytes <- function(key) {
        if (!length(io)) return(0)
        values <- unlist(strsplit(io, " ", fixed = TRUE), use.names = FALSE)
        values <- values[startsWith(values, paste0(key, "="))]
        sum(as.numeric(sub(paste0("^", key, "="), "", values)))
    }
    if (length(usage) != 1L || !is.finite(usage)) cli::cli_abort("CPU accounting unavailable")
    list(cpu_seconds = usage / 1e6, read_bytes = bytes("rbytes"), written_bytes = bytes("wbytes"))
}

#' Verify that a known worker subtree has no remaining tasks after stop.
#' @param application Managed application.
#' @param unit Owned service beneath that application.
#' @return Whether the kernel subtree is gone or empty.
#' @noRd
workflowServiceEmpty <- function(application, unit) {
    path <- file.path(application$cgroup, unit)
    if (!dir.exists(path)) return(TRUE)
    files <- c(file.path(path, "cgroup.procs"), list.files(path, "^cgroup[.]procs$",
        recursive = TRUE, full.names = TRUE))
    all(vapply(unique(files), function(file) file.exists(file) &&
        !length(readLines(file, warn = FALSE)), logical(1)))
}

#' Identify one live controller across PID reuse and host restarts.
#' @return Small process identity; non-Linux legacy callers have no managed identity.
#' @noRd
workflowExecutionProcess <- function(pid = Sys.getpid()) {
    if (!is.numeric(pid) || length(pid) != 1L || !is.finite(pid) ||
        pid < 1 || pid != floor(pid)) stop("invalid process identity")
    if (!identical(Sys.info()[["sysname"]], "Linux")) {
        return(list(pid = as.integer(pid), platform = Sys.info()[["sysname"]]))
    }
    path <- file.path("/proc", as.character(as.integer(pid)), "stat")
    stat <- strsplit(sub("^.*[)] ", "", readLines(path, warn = FALSE)), " ",
        fixed = TRUE)[[1L]]
    list(pid = as.integer(pid), start = stat[[20L]],
        boot = readLines("/proc/sys/kernel/random/boot_id", warn = FALSE))
}

#' Check project identity for cleanup without writing into a replacement project.
#' @param store Original owner.
#' @return Whether the original store is still current.
#' @noRd
workflowExecutionStoreValid <- function(store) {
    tryCatch({ checkWorkflowEntityStore(store); TRUE }, saved_entity_scope = function(error) FALSE,
        saved_entity_invalid = function(error) FALSE, saved_entity_corrupt = function(error) FALSE)
}

#' Run one fixed diagnostic controller-death case inside the managed application.
#' @param application Diagnostic resource owner; never exposed as a science adapter.
#' @return No normal result when the test controller is killed at its marker.
#' @noRd
workflowControllerFaultMain <- function(application) {
    config <- application$config
    if (!isTRUE(config$policy$diagnostic) || !config$fault_stage %in%
        c("before_head_publication", "after_head_publication")) stop("invalid diagnostic fault case")
    project <- normalizePath(config$fault_project, mustWork = TRUE)
    if (!startsWith(project, paste0(config$root, "/"))) stop("fault project escapes owned application")
    identity <- readSavedStoreIdentity(file.path(project, ".workflow-entities", "store.json"))
    store <- newWorkflowEntityStore(project, identity$project_id, identity$workflow_id)
    manager <- newDIAStateAdapter(store, reopen = TRUE)
    owner <- newWorkflowExecutionController(manager, application)
    on.exit(closeDIAStateAdapter(manager))
    testthat::local_mocked_bindings(diaStateBoundary = function(stage) {
        if (identical(stage, config$fault_stage)) {
            writeWorkflowExecutionMessage(config$fault_marker, list(stage = stage,
                request_id = owner$job$request$token$request_id))
            Sys.sleep(60)
            stop("fault controller was not killed by the test")
        }
    }, .package = "MultiScholaR")
    runWorkflowAction(workflowFixtureOperation, list(state = getDIAStateReference(manager)), list(mode = "identity"), owner)
    repeat {
        value <- pollWorkflowAction(owner)
        if (!value$status %in% c("admitting", "running")) stop("fault boundary was not reached")
        Sys.sleep(0.02)
    }
}

#' Create the accounted application parent before starting any controller.
#' @param stem Internally generated application unit stem.
#' @param policy Fixed server-owned resource policy.
#' @return Actual cgroup path reported by the user manager.
#' @noRd
createWorkflowApplicationSlice <- function(stem, policy) {
    if (!grepl("^multischolarwfx[a-f0-9]{16}$", stem) ||
        !identical(policy, workflowExecutionPolicy(policy$diagnostic))) stop("invalid application slice")
    properties <- c("Description", "s", "MultiScholaR guarded application", "MemoryAccounting", "b", "true",
        "CPUAccounting", "b", "true", "IOAccounting", "b", "true",
        "MemoryMax", "t", sprintf("%.0f", policy$memory_max), "MemorySwapMax", "t", "0")
    slice <- paste0(stem, ".slice")
    workflowSystemd(c("call", "org.freedesktop.systemd1", "/org/freedesktop/systemd1",
        "org.freedesktop.systemd1.Manager", "StartTransientUnit", "ssa(sv)a(sa(sv))",
        slice, "fail", "6", properties, "0"), command = "busctl")
    workflowServiceProperties(slice)$ControlGroup
}


#' Validate parameters for the finite internal operation registry.
#' @param operation Registered operation identifier.
#' @param parameters Bounded immutable parameters.
#' @param policy Server-owned resource policy.
#' @param resolved Whether native precedence has already produced complete parameters.
#' @return NULL invisibly or explicit refusal before dispatch.
#' @noRd
validateWorkflowActionParameters <- function(operation, parameters, policy, resolved = TRUE) {
    entry <- workflowActionSpec(operation)
    if (!is.list(parameters) || hasMutableWorkflowEntityOwner(parameters) ||
        as.numeric(object.size(parameters)) > entry$parameter_limit || anyDuplicated(names(parameters)))
        cli::cli_abort("unsupported execution parameters", class = "workflow_execution_invalid")
    # The entry names its validator; it is resolved here so a sourced copy of this file dispatches its own.
    get(entry$validate, mode = "function")(parameters, policy, resolved)
    invisible(NULL)
}

#' P14 takes no caller-selectable settings; the native artifact must be present.
#' @noRd
validateWorkflowPrecursorRollupParameters <- function(parameters, policy, resolved = TRUE) {
    if (length(parameters)) cli::cli_abort("P14 has no caller-selectable execution settings")
    diaRollupNativeInfo()
    invisible(NULL)
}

#' The sample filter's cutoff and inclusion list, complete once resolved.
#' @noRd
validateWorkflowSampleFilterParameters <- function(parameters, policy, resolved = TRUE) {
    allowed <- c("peptides_per_sample_cutoff", "inclusion_list")
    if (length(parameters) && (is.null(names(parameters)) ||
        anyNA(names(parameters)) || any(!names(parameters) %in% allowed)) ||
        resolved && !setequal(names(parameters), allowed))
        cli::cli_abort("sample filter parameters required or unknown parameter supplied")
    k <- parameters$peptides_per_sample_cutoff
    if ((resolved || !is.null(k)) && (!is.numeric(k) || length(k) != 1L ||
        !is.finite(k) || k < 1 || k > .Machine$integer.max || k != floor(k)))
        cli::cli_abort("cutoff must be a positive integer")
    if (!is.null(parameters$inclusion_list) && !is.character(parameters$inclusion_list))
        cli::cli_abort("inclusion entries must be character")
    invisible(NULL)
}

#' The imputation's proportion threshold and optional design exclusion column (KIT-009).
#' @noRd
validateWorkflowPeptideImputationParameters <- function(parameters, policy, resolved = TRUE) {
    allowed <- c("proportion_missing_values", "exclusion_column", "skip")
    if (length(parameters) && (is.null(names(parameters)) || anyNA(names(parameters)) ||
        any(!names(parameters) %in% allowed)) || resolved && !"proportion_missing_values" %in% names(parameters))
        cli::cli_abort("imputation parameters required or unknown parameter supplied")
    # KIT-049: the declared skip is a logical scalar.
    if (!is.null(parameters$skip) && (!is.logical(parameters$skip) || length(parameters$skip) != 1L || is.na(parameters$skip)))
        cli::cli_abort("skip must be TRUE or FALSE")
    p <- parameters$proportion_missing_values
    if ((resolved || !is.null(p)) && (!is.numeric(p) || length(p) != 1L || !is.finite(p) || p < 0 || p > 1))
        cli::cli_abort("proportion of missing values must lie in [0, 1]")
    exclusion <- parameters$exclusion_column
    if (!is.null(exclusion) && (!is.character(exclusion) || length(exclusion) != 1L))
        cli::cli_abort("exclusion column must be a single name")
    invisible(NULL)
}

#' The Q-value filter's fixed confidence contract (KIT-014): three thresholds in (0, 1], the proteotypic
#' flag and the provenance mode.
#' @noRd
validateWorkflowPeptideQvalueParameters <- function(parameters, policy, resolved = TRUE) {
    allowed <- c("qvalue_threshold", "global_qvalue_threshold", "global_pg_qvalue_threshold", "proteotypic_only",
        "confidence_provenance_mode")
    if (length(parameters) && (is.null(names(parameters)) || anyNA(names(parameters)) ||
        any(!names(parameters) %in% allowed)) || resolved && !setequal(names(parameters), allowed))
        cli::cli_abort("Q-value parameters required or unknown parameter supplied")
    for (name in allowed[1:3]) {
        value <- parameters[[name]]
        if ((resolved || !is.null(value)) && (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
            value <= 0 || value > 1))
            cli::cli_abort("confidence thresholds must lie in (0, 1]")
    }
    flag <- parameters$proteotypic_only
    if ((resolved || !is.null(flag)) && (!is.logical(flag) || length(flag) != 1L || is.na(flag)))
        cli::cli_abort("proteotypic flag must be TRUE or FALSE")
    mode <- parameters$confidence_provenance_mode
    if ((resolved || !is.null(mode)) && (!is.character(mode) || length(mode) != 1L ||
        !mode %in% c("diann_main_report", "upstream_prefiltered")))
        cli::cli_abort("unsupported confidence provenance mode")
    invisible(NULL)
}

#' The protein identification evidence filter's two positive integer thresholds (KIT-014).
#' @noRd
validateWorkflowPeptideProteinEvidenceParameters <- function(parameters, policy, resolved = TRUE) {
    allowed <- c("num_peptides_per_protein_thresh", "num_peptidoforms_per_protein_thresh")
    if (length(parameters) && (is.null(names(parameters)) || anyNA(names(parameters)) ||
        any(!names(parameters) %in% allowed)) || resolved && !setequal(names(parameters), allowed))
        cli::cli_abort("protein evidence parameters required or unknown parameter supplied")
    for (name in allowed) {
        value <- parameters[[name]]
        if ((resolved || !is.null(value)) && (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
            value < 1 || value != floor(value)))
            cli::cli_abort("evidence thresholds must be positive integers")
    }
    invisible(NULL)
}

#' The peptide intensity filter's strict flag, replicate and group minimums, percentile and optional
#' grouping variable (KIT-014).
#' @noRd
validateWorkflowPeptideIntensityParameters <- function(parameters, policy, resolved = TRUE) {
    allowed <- c("strict_mode", "min_reps_per_group", "min_groups", "peptides_intensity_cutoff_percentile", "grouping_variable")
    required <- allowed[1:4]
    if (length(parameters) && (is.null(names(parameters)) || anyNA(names(parameters)) ||
        any(!names(parameters) %in% allowed)) || resolved && !all(required %in% names(parameters)))
        cli::cli_abort("intensity filter parameters required or unknown parameter supplied")
    strict <- parameters$strict_mode
    if ((resolved || !is.null(strict)) && (!is.logical(strict) || length(strict) != 1L || is.na(strict)))
        cli::cli_abort("strict mode must be TRUE or FALSE")
    for (name in c("min_reps_per_group", "min_groups")) {
        value <- parameters[[name]]
        if ((resolved || !is.null(value)) && (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
            value < 1 || value != floor(value)))
            cli::cli_abort("replicate and group minimums must be positive integers")
    }
    percentile <- parameters$peptides_intensity_cutoff_percentile
    if ((resolved || !is.null(percentile)) && (!is.numeric(percentile) || length(percentile) != 1L ||
        !is.finite(percentile) || percentile < 0 || percentile > 100))
        cli::cli_abort("intensity percentile must lie in [0, 100]")
    grouping <- parameters$grouping_variable
    if (!is.null(grouping) && (!is.character(grouping) || length(grouping) != 1L || is.na(grouping) || !nzchar(grouping)))
        cli::cli_abort("grouping variable must be a single column name")
    invisible(NULL)
}

#' The peptide replicate filter's replicate group column (KIT-014).
#' @noRd
validateWorkflowPeptideReplicateParameters <- function(parameters, policy, resolved = TRUE) {
    allowed <- "replicate_group_column"
    if (length(parameters) && (is.null(names(parameters)) || anyNA(names(parameters)) ||
        any(!names(parameters) %in% allowed)) || resolved && !setequal(names(parameters), allowed))
        cli::cli_abort("replicate filter parameters required or unknown parameter supplied")
    column <- parameters$replicate_group_column
    if ((resolved || !is.null(column)) && (!is.character(column) || length(column) != 1L || is.na(column) || !nzchar(column)))
        cli::cli_abort("replicate group column must be a single column name")
    invisible(NULL)
}

#' The IQ rollup declares no scientific parameters; the working directories travel outside the identity
#' (KIT-015).
#' @noRd
validateWorkflowProteinIQRollupParameters <- function(parameters, policy, resolved = TRUE) {
    allowed <- c("peptide_qc_dir", "protein_qc_dir")
    if (length(parameters) && (is.null(names(parameters)) || anyNA(names(parameters)) || any(!names(parameters) %in% allowed)))
        cli::cli_abort("unknown IQ rollup parameter supplied")
    for (name in allowed) {
        value <- parameters[[name]]
        if (!is.null(value) && (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)))
            cli::cli_abort("IQ working directories must be single paths")
    }
    invisible(NULL)
}

#' The accession cleanup's delimiter, aggregation method and optional FASTA table path (KIT-015).
#' @noRd
validateWorkflowProteinAccessionParameters <- function(parameters, policy, resolved = TRUE) {
    allowed <- c("delimiter", "aggregation_method", "fasta_table", "fasta_sha256")
    required <- allowed[1:2]
    if (length(parameters) && (is.null(names(parameters)) || anyNA(names(parameters)) ||
        any(!names(parameters) %in% allowed)) || resolved && !all(required %in% names(parameters)))
        cli::cli_abort("accession cleanup parameters required or unknown parameter supplied")
    for (name in allowed) {
        value <- parameters[[name]]
        if ((resolved && name %in% required || !is.null(value)) &&
            (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)))
            cli::cli_abort("accession cleanup parameters must be single strings")
    }
    invisible(NULL)
}

#' The DIA_LIMPA peptide stage's parameters (KIT-047): the log2 matrix declares none; the peptide accession cleanup
#' takes the protein cleanup's strings plus the FASTA accession column; the DPC quantification its engine, slope and chunk.
#' @noRd
validateWorkflowPeptideLog2Parameters <- function(parameters, policy, resolved = TRUE) {
    if (length(parameters)) cli::cli_abort("the log2 matrix step takes no parameters")
    invisible(NULL)
}

#' @noRd
validateWorkflowPeptideAccessionParameters <- function(parameters, policy, resolved = TRUE) {
    allowed <- c("delimiter", "aggregation_method", "seqinr_accession_column", "fasta_table", "fasta_sha256")
    required <- allowed[1:2]
    if (length(parameters) && (is.null(names(parameters)) || anyNA(names(parameters)) ||
        any(!names(parameters) %in% allowed)) || resolved && !all(required %in% names(parameters)))
        cli::cli_abort("peptide accession cleanup parameters required or unknown parameter supplied")
    for (name in allowed) {
        value <- parameters[[name]]
        if ((resolved && name %in% required || !is.null(value)) &&
            (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)))
            cli::cli_abort("peptide accession cleanup parameters must be single strings")
    }
    if (!is.null(parameters$aggregation_method) && !parameters$aggregation_method %in% c("sum", "mean", "median"))
        cli::cli_abort("aggregation method must be sum, mean or median")
    invisible(NULL)
}

#' @noRd
validateWorkflowProteinDpcQuantParameters <- function(parameters, policy, resolved = TRUE) {
    allowed <- c("engine", "dpc_slope", "chunk")
    if (length(parameters) && (is.null(names(parameters)) || anyNA(names(parameters)) ||
        any(!names(parameters) %in% allowed)) || resolved && !setequal(names(parameters), allowed))
        cli::cli_abort("DPC quantification parameters required or unknown parameter supplied")
    engine <- parameters$engine
    if ((resolved || !is.null(engine)) && (!is.character(engine) || length(engine) != 1L || !engine %in% c("fastlimpadpc", "limpa")))
        cli::cli_abort("DPC engine must be fastlimpadpc or limpa")
    slope <- parameters$dpc_slope
    if ((resolved || !is.null(slope)) && (!is.numeric(slope) || length(slope) != 1L || !is.finite(slope) || slope <= 0))
        cli::cli_abort("DPC slope must be a positive number")
    chunk <- parameters$chunk
    if ((resolved || !is.null(chunk)) && (!is.numeric(chunk) || length(chunk) != 1L || !is.finite(chunk) || chunk < 1 || chunk != trunc(chunk)))
        cli::cli_abort("DPC chunk must be a positive whole number")
    invisible(NULL)
}

#' The protein intensity filter's strict flag, replicate and group minimums and percentile (KIT-015).
#' @noRd
validateWorkflowProteinIntensityParameters <- function(parameters, policy, resolved = TRUE) {
    allowed <- c("strict_mode", "min_reps_per_group", "min_groups", "proteins_intensity_cutoff_percentile")
    if (length(parameters) && (is.null(names(parameters)) || anyNA(names(parameters)) ||
        any(!names(parameters) %in% allowed)) || resolved && !setequal(names(parameters), allowed))
        cli::cli_abort("protein intensity filter parameters required or unknown parameter supplied")
    strict <- parameters$strict_mode
    if ((resolved || !is.null(strict)) && (!is.logical(strict) || length(strict) != 1L || is.na(strict)))
        cli::cli_abort("strict mode must be TRUE or FALSE")
    for (name in c("min_reps_per_group", "min_groups")) {
        value <- parameters[[name]]
        if ((resolved || !is.null(value)) && (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
            value < 1 || value != floor(value)))
            cli::cli_abort("replicate and group minimums must be positive integers")
    }
    percentile <- parameters$proteins_intensity_cutoff_percentile
    if ((resolved || !is.null(percentile)) && (!is.numeric(percentile) || length(percentile) != 1L ||
        !is.finite(percentile) || percentile < 0 || percentile > 100))
        cli::cli_abort("intensity percentile must lie in [0, 100]")
    invisible(NULL)
}

#' The duplicate removal's aggregation function name (KIT-015).
#' @noRd
validateWorkflowProteinDuplicateParameters <- function(parameters, policy, resolved = TRUE) {
    allowed <- "aggregation_method"
    if (length(parameters) && (is.null(names(parameters)) || anyNA(names(parameters)) ||
        any(!names(parameters) %in% allowed)) || resolved && !setequal(names(parameters), allowed))
        cli::cli_abort("duplicate removal parameters required or unknown parameter supplied")
    method <- parameters$aggregation_method
    if ((resolved || !is.null(method)) && (!is.character(method) || length(method) != 1L || is.na(method) ||
        !method %in% c("mean", "median", "sum", "max", "min")))
        cli::cli_abort("aggregation method must be one of mean, median, sum, max, min")
    invisible(NULL)
}

#' The protein replicate filter's grouping variable and the recorded core count (KIT-010).
#' @noRd
validateWorkflowProteinReplicateParameters <- function(parameters, policy, resolved = TRUE) {
    allowed <- c("grouping_variable", "parallel_cores")
    if (length(parameters) && (is.null(names(parameters)) || anyNA(names(parameters)) ||
        any(!names(parameters) %in% allowed)) || resolved && !"grouping_variable" %in% names(parameters))
        cli::cli_abort("replicate filter parameters required or unknown parameter supplied")
    g <- parameters$grouping_variable
    if ((resolved || !is.null(g)) && (!is.character(g) || length(g) != 1L || is.na(g) || !nzchar(g)))
        cli::cli_abort("grouping variable must be a single column name")
    cores <- parameters$parallel_cores
    if (!is.null(cores) && (!is.numeric(cores) || length(cores) != 1L || !is.finite(cores) || cores < 1 || cores != floor(cores)))
        cli::cli_abort("parallel cores must be a positive integer")
    invisible(NULL)
}

#' The diagnostic fixture's modes, under a diagnostic policy only.
#' @noRd
validateWorkflowFixtureParameters <- function(parameters, policy, resolved = TRUE) {
    if (!isTRUE(policy$diagnostic)) cli::cli_abort("fixture requires diagnostic policy")
    if (is.null(names(parameters)) || any(!names(parameters) %in% c("mode", "seconds", "bytes")) ||
        !is.character(parameters$mode) || length(parameters$mode) != 1L || is.na(parameters$mode) ||
        !parameters$mode %in% c("identity", "fail", "sleep", "allocate", "forbidden_commit"))
        cli::cli_abort("unsupported execution parameters", class = "workflow_execution_invalid")
    if (parameters$mode == "sleep" && (!is.numeric(parameters$seconds) || length(parameters$seconds) != 1L ||
        !is.finite(parameters$seconds) || parameters$seconds < 0 || parameters$seconds > 60))
        cli::cli_abort("invalid fixture delay")
    if (parameters$mode == "allocate" && (!is.numeric(parameters$bytes) || length(parameters$bytes) != 1L ||
        !is.finite(parameters$bytes) || parameters$bytes < 8 || parameters$bytes > 2 * 1024^3 || parameters$bytes %% 8 != 0))
        cli::cli_abort("invalid bounded allocation fixture")
    invisible(NULL)
}

#' Compute one frozen scientific checkpoint inside W, without a mutation lease.
#' @param request Immutable action request.
#' @return Existing bounded result/checkpoint structures.
#' @noRd
computeWorkflowPeptideAction <- function(request) {
    store <- newWorkflowEntityStore(request$project_root, request$owner$project_id, request$owner$workflow_id)
    on.exit(closeWorkflowEntityStore(store))
    if (!identical(store$identity, request$owner)) cli::cli_abort("project incarnation changed")
    manager <- newDIAStateAdapter(store, reopen = TRUE)
    on.exit(closeDIAStateAdapter(manager), add = TRUE, after = FALSE)
    validateWorkflowActionParameters(request$operation, request$parameters, request$config$policy)
    validateWorkflowActionInput(request$operation, manager, request$input_references)
    variant <- workflowQualificationVariant(request$config, request$operation)
    if (!is.null(request$config$qualification) && !identical(request$token$qualification,
        digest::digest(request$config$qualification, algo = "sha256")))
        cli::cli_abort("qualification variant token differs", class = "workflow_execution_invalid")
    entry <- workflowActionSpec(request$operation)
    if (!entry$compute %in% c("peptide", "protein")) cli::cli_abort("scientific adapter unavailable")
    value <- get(entry$handler, mode = "function")(manager, request, variant)
    value$operation <- request$operation
    value
}

#' Compute handlers named by the registry: each runs the unchanged reference entry point of its
#' declaration with the dispatched worker cap and the qualification variant.
#' @noRd
computeWorkflowPrecursorRollupAction <- function(manager, request, variant) {
    computeDIAReferencePrecursorRollup(manager, request$input_references$state, worker_cap = request$worker_cap)
}

#' @noRd
computeWorkflowSampleFilterAction <- function(manager, request, variant) {
    computeDIAReferenceSampleFilter(manager, request$input_references$state, request$parameters,
        worker_cap = request$worker_cap, output_mode = variant)
}

#' The experiment paths of a request's publication context, when the controller wrote one beside the
#' request (a GUI action), else NULL (KIT-010: the worker records the project's output path in the state
#' configuration as the GUI does).
#' @noRd
workflowRequestExperimentPaths <- function(request) {
    path <- if (is.character(request$directory) && length(request$directory) == 1L)
        file.path(request$directory, "publication-context.rds") else NULL
    if (is.null(path) || !file.exists(path)) return(NULL)
    publication <- tryCatch(readRDS(path), error = function(error) NULL)
    if (!is.list(publication)) return(NULL)
    if (is.list(publication$experiment_paths)) return(publication$experiment_paths)
    if (is.list(publication$proteomics)) return(publication$proteomics)
    NULL
}

#' @noRd
computeWorkflowProteinReplicateAction <- function(manager, request, variant) {
    paths <- workflowRequestExperimentPaths(request)
    if (!is.null(paths)) manager$.context$experiment_paths <- paths
    computeDIAReferenceProteinReplicateFilter(manager, request$input_references$state, request$parameters,
        worker_cap = request$worker_cap, variant = variant)
}

#' The KIT-015 protein QC handlers: the worker's request directory holds IQ's working files until the
#' controller's publication phase copies them; the accession cleanup reads the FASTA table from the bound
#' path or the owner's source directory.
#' @noRd
computeWorkflowProteinIQRollupAction <- function(manager, request, variant) {
    paths <- workflowRequestExperimentPaths(request)
    if (!is.null(paths)) manager$.context$experiment_paths <- paths
    parameters <- request$parameters
    if (is.null(parameters$peptide_qc_dir)) parameters$peptide_qc_dir <- file.path(request$directory, "iq", "peptide")
    if (is.null(parameters$protein_qc_dir)) parameters$protein_qc_dir <- file.path(request$directory, "iq", "protein")
    computeDIAReferenceIQRollup(manager, request$input_references$state, parameters,
        worker_cap = request$worker_cap, variant = variant)
}

#' @noRd
computeWorkflowProteinAccessionAction <- function(manager, request, variant) {
    paths <- workflowRequestExperimentPaths(request)
    if (!is.null(paths)) manager$.context$experiment_paths <- paths
    computeDIAReferenceAccessionCleanup(manager, request$input_references$state, request$parameters,
        worker_cap = request$worker_cap, variant = variant)
}

#' KIT-047: the DIA_LIMPA peptide stage's handlers; the accession cleanup reads the FASTA table as the protein one does.
#' @noRd
computeWorkflowPeptideLog2Action <- function(manager, request, variant) {
    computeDIAReferencePeptideLog2Matrix(manager, request$input_references$state, request$parameters,
        worker_cap = request$worker_cap, variant = variant)
}

#' @noRd
computeWorkflowPeptideAccessionAction <- function(manager, request, variant) {
    paths <- workflowRequestExperimentPaths(request)
    if (!is.null(paths)) manager$.context$experiment_paths <- paths
    computeDIAReferencePeptideAccessionCleanup(manager, request$input_references$state, request$parameters,
        worker_cap = request$worker_cap, variant = variant)
}

#' @noRd
computeWorkflowProteinDpcQuantAction <- function(manager, request, variant) {
    paths <- workflowRequestExperimentPaths(request)
    if (!is.null(paths)) manager$.context$experiment_paths <- paths
    computeDIAReferenceProteinDpcQuant(manager, request$input_references$state, request$parameters,
        worker_cap = request$worker_cap, variant = variant)
}

#' KIT-041: the DIA_LIMPA protein stage's handlers and validators.
#' @noRd
computeWorkflowProteinNormalisationAction <- function(manager, request, variant) {
    computeDIAReferenceProteinNormalisation(manager, request$input_references$state, request$parameters,
        worker_cap = request$worker_cap, variant = variant)
}

#' @noRd
computeWorkflowProteinPeptideCountAction <- function(manager, request, variant) {
    computeDIAReferenceProteinPeptideCountFilter(manager, request$input_references$state, request$parameters,
        worker_cap = request$worker_cap, variant = variant)
}

#' @noRd
validateWorkflowProteinNormalisationParameters <- function(parameters, policy, resolved = TRUE) {
    allowed <- "norm_method"
    if (length(parameters) && (is.null(names(parameters)) || anyNA(names(parameters)) ||
        any(!names(parameters) %in% allowed)) || resolved && !setequal(names(parameters), allowed))
        cli::cli_abort("protein normalisation parameters required or unknown parameter supplied")
    method <- parameters$norm_method
    if ((resolved || !is.null(method)) && (!is.character(method) || length(method) != 1L ||
        !method %in% c("none", "cyclicloess", "quantile", "scale")))
        cli::cli_abort("normalisation method must be none, cyclicloess, quantile or scale")
    invisible(NULL)
}

#' @noRd
validateWorkflowProteinPeptideCountParameters <- function(parameters, policy, resolved = TRUE) {
    allowed <- c("num_peptides_per_protein_thresh", "num_peptidoforms_per_protein_thresh")
    if (length(parameters) && (is.null(names(parameters)) || anyNA(names(parameters)) ||
        any(!names(parameters) %in% allowed)) || resolved && !setequal(names(parameters), allowed))
        cli::cli_abort("peptide count filter parameters required or unknown parameter supplied")
    for (name in allowed) {
        value <- parameters[[name]]
        if ((resolved || !is.null(value)) && (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value < 0 || value != trunc(value)))
            cli::cli_abort("peptide count thresholds must be non-negative whole numbers")
    }
    invisible(NULL)
}

#' @noRd
computeWorkflowProteinIntensityAction <- function(manager, request, variant) {
    computeDIAReferenceProteinIntensityFilter(manager, request$input_references$state, request$parameters,
        worker_cap = request$worker_cap, variant = variant)
}

#' @noRd
computeWorkflowProteinDuplicateAction <- function(manager, request, variant) {
    computeDIAReferenceProteinDuplicateRemoval(manager, request$input_references$state, request$parameters,
        worker_cap = request$worker_cap, variant = variant)
}

#' An in-process declaration's parameters are validated by its own runner, never by a worker request (KIT-024).
#' @noRd
validateWorkflowInProcessParameters <- function(parameters, policy, resolved = TRUE) invisible(NULL)

#' The DA publication's one parameter (KIT-019): the annotation table as a saved reference, or absent.
#' @noRd
validateWorkflowDaPublicationParameters <- function(parameters, policy, resolved = TRUE) {
    if (length(parameters) && (is.null(names(parameters)) || anyNA(names(parameters)) ||
        any(!names(parameters) %in% "annotations")))
        cli::cli_abort("DA publication takes only an annotations reference")
    annotations <- parameters$annotations
    if (!is.null(annotations) && !isSavedWorkflowReference(annotations))
        cli::cli_abort("annotations must be a saved table reference")
    invisible(NULL)
}

#' The DA publication in the worker (KIT-019): the accepted manifest from the head snapshot, the annotation
#' table from its reference, the files staged under the request directory (`da-publication/<role>`); the
#' publication kernels resolve their tables through the active owner, which the state scope names (KIT-022).
#' The pages are not saved as assets here (the controller holds the head).
#' @return The completion value: the head checkpoint and reference (nothing new is committed) and the inventory.
#' @noRd
computeWorkflowDaPublicationAction <- function(manager, request, variant) {
    paths <- workflowRequestExperimentPaths(request)
    if (!is.null(paths)) manager$.context$experiment_paths <- paths
    head <- getDIAStateReference(manager)
    if (!identical(head$entity_id, request$input_references$state$entity_id))
        cli::cli_abort("DA publication input is not the head state", class = "workflow_execution_stale")
    manifest_reference <- manager$.context$snapshot$da_results
    if (is.null(manifest_reference)) cli::cli_abort("no accepted DA manifest", class = "workflow_execution_invalid")
    manifest <- resolveDIAStateValue(manager$.context, manifest_reference)
    # KIT-022: the manifest's source is the head (the DA step commits its arguments before the contrasts).
    head_record <- readDIAStateMetadata(manager, head)$data
    if (!isDIAResultManifest(manifest) || !identical(manifest$source$entity_id, head$entity_id))
        cli::cli_abort("DA manifest does not belong to the head state", class = "workflow_execution_stale")
    # As the in-process path: no annotation table means an empty one.
    annotations <- if (!is.null(request$parameters$annotations))
        materializeDIATable(manager$.context, request$parameters$annotations) else data.frame()
    root <- file.path(request$directory, "da-publication")
    staged <- stats::setNames(lapply(c("da_output_dir", "publication_graphs_dir"), function(role) {
        path <- file.path(root, role)
        dir.create(path, recursive = TRUE, showWarnings = FALSE)
        path
    }), c("da_output_dir", "publication_graphs_dir"))
    # KIT-022: the state scope names its owner for the publication kernels' readers; no session is needed.
    withDIAStateScope(manager, function()
        produceDIADaPublication(manager, manifest, annotations, manifest$parameters$qvalue, staged, save_assets = FALSE))
    files <- list.files(root, recursive = TRUE)
    list(checkpoint = head_record, reference = head, state_reference = head,
        state_name = manager$current_state, manifest_reference = manifest_reference,
        summary = list(files = length(files), bytes = sum(file.size(file.path(root, files)))), files = files)
}

#' Complete a DA publication (KIT-019): nothing new is committed (the manifest was committed before the request);
#' the head is republished unchanged so the controller's acceptance records the head state as the result.
#' @noRd
completeDIADaPublication <- function(workflow, value) {
    manager <- workflow$state_manager
    head <- getDIAStateReference(manager)
    if (!identical(head$entity_id, value$reference$entity_id))
        cli::cli_abort("DA publication result belongs to a superseded head", class = "workflow_execution_stale")
    context <- manager$.context
    publishDIAStateHead(context, context$snapshot)
    list(stateReference = head, manifest_reference = value$manifest_reference)
}

#' Publish the worker's staged DA files into the project (KIT-019) through the owner's action-output
#' journal, and save each contrast's page as the KIT-016 rendered asset (the controller holds the head).
#' @noRd
publishWorkflowDaResult <- function(owner, paths, result) {
    manager <- owner$manager
    entry <- workflowActionSpec(result$operation)
    reference <- result$state_reference
    context_path <- file.path(paths$request_dir, "publication-context.rds")
    publication <- if (file.exists(context_path)) readRDS(context_path) else NULL
    experiment_paths <- if (is.list(publication) && is.list(publication$experiment_paths)) publication$experiment_paths else
        if (is.list(publication) && is.list(publication$proteomics)) publication$proteomics else {
            root <- file.path(paths$request_dir, "required-output")
            dir.create(root, showWarnings = FALSE)
            setupDirectories(root, entry$omic_type)$proteomics
        }
    if (!is.list(experiment_paths) || is.null(experiment_paths$da_output_dir) || is.null(experiment_paths$publication_graphs_dir))
        cli::cli_abort("publication context is invalid", class = "workflow_execution_invalid")
    manager$.context$experiment_paths <- experiment_paths
    manifest_reference <- manager$.context$snapshot$da_results
    if (is.null(manifest_reference)) cli::cli_abort("no accepted DA manifest", class = "workflow_execution_invalid")
    manifest <- resolveDIAStateValue(manager$.context, manifest_reference)
    staged_root <- file.path(paths$request_dir, "da-publication")
    workflowExecutionBoundary("publication")
    publishDIAActionOutputs(manager, reference, "da", experiment_paths[c("da_output_dir", "publication_graphs_dir")],
        function(staged) {
            for (role in names(staged)) {
                members <- list.files(file.path(staged_root, role), full.names = TRUE, all.files = TRUE, no.. = TRUE)
                if (length(members) && !all(file.copy(members, staged[[role]], recursive = TRUE)))
                    cli::cli_abort("cannot stage the worker's DA files", class = "saved_entity_outputs_pending")
            }
            pages <- file.path(staged$publication_graphs_dir, "Interactive_Volcano_Plots",
                paste0(gsub("[^A-Za-z0-9_-]", "_", manifest$contrasts$friendly_names), ".html"))
            assertDIADaPublicationFiles(manifest, staged, pages)
            for (i in seq_along(pages)) saveDIAInteractiveVolcanoAsset(manager, manifest, manifest$contrasts$friendly_names[[i]],
                manifest$parameters$qvalue, readBin(pages[[i]], "raw", n = file.size(pages[[i]])))
            list(completed = TRUE)
        }, identity = diaDaPublicationIdentity(manifest))
    invisible(NULL)
}

#' Recover the accepted DA manifest for the module after the managed publication (KIT-019).
#' @noRd
readManagedDaResult <- function(manager, result, workflow_data) {
    manifest_reference <- manager$.context$snapshot$da_results
    if (is.null(manifest_reference)) cli::cli_abort("no accepted DA manifest", class = "workflow_execution_invalid")
    manifest <- resolveDIAStateValue(manager$.context, manifest_reference)
    list(stateReference = result$state_reference, manifest = manifest,
        resultText = diaDaStatusText(nrow(manifest$contrasts)))
}

#' @noRd
computeWorkflowPeptideImputationAction <- function(manager, request, variant) {
    computeDIAReferencePeptideGroups(manager, request$input_references$state, "peptideMissingValueImputation",
        request$parameters, worker_cap = request$worker_cap, variant = variant)
}

#' The KIT-014 peptide QC handlers: each section of the peptide-group family through its compute-only entry
#' point, and the Q-value filter through its own.
#' @noRd
computeWorkflowPeptideProteinEvidenceAction <- function(manager, request, variant) {
    computeDIAReferencePeptideGroups(manager, request$input_references$state, "filterMinNumPeptidesPerProtein",
        request$parameters, worker_cap = request$worker_cap, variant = variant)
}

#' @noRd
computeWorkflowPeptideIntensityAction <- function(manager, request, variant) {
    computeDIAReferencePeptideGroups(manager, request$input_references$state, "peptideIntensityFiltering",
        request$parameters, worker_cap = request$worker_cap, variant = variant)
}

#' @noRd
computeWorkflowPeptideReplicateAction <- function(manager, request, variant) {
    computeDIAReferencePeptideGroups(manager, request$input_references$state, "removePeptidesWithOnlyOneReplicate",
        request$parameters, worker_cap = request$worker_cap, variant = variant)
}

#' @noRd
computeWorkflowPeptideQvalueAction <- function(manager, request, variant) {
    computeDIAReferenceQvalue(manager, request$input_references$state, request$parameters,
        worker_cap = request$worker_cap, variant = variant)
}

#' Resolve only server-owned finite qualification variants for installed actions.
#' @param config Frozen application configuration, never GUI scientific parameters.
#' @param operation Existing registered action identifier.
#' @return Legacy default or an explicitly bound experimental sample-filter mode.
#' @noRd
workflowQualificationVariant <- function(config, operation) {
    qualification <- config$qualification
    # KIT-008: an unqualified launch (production, managed or not) runs the registry's production variant.
    if (is.null(qualification)) return(workflowActionSpec(operation)$default_variant %||% "legacy")
    if (!identical(config$mode, "adapters") ||
        !is.list(qualification) || !identical(names(qualification),
            c("programme", "variant", "member_sha256", "deadline")) ||
        !identical(qualification$programme, "perf106-repair-v2") ||
        !is.character(qualification$member_sha256) || length(qualification$member_sha256) != 1L ||
        is.na(qualification$member_sha256) || !grepl("^[a-f0-9]{64}$", qualification$member_sha256) ||
        !is.numeric(qualification$deadline) || length(qualification$deadline) != 1L ||
        is.na(qualification$deadline) || !qualification$deadline %in% c(300, 600, 1800) ||
        !is.character(qualification$variant) || length(qualification$variant) != 1L ||
        is.na(qualification$variant) || !qualification$variant %in% workflowActionSpec(operation)$variants)
        cli::cli_abort("unsupported qualification variant", class = "workflow_execution_invalid")
    qualification$variant
}

#' Publish accepted peptide QC figures using existing checkpoint/output recovery.
#' @param owner Leased controller.
#' @param paths Request-owned publication paths.
#' @param result Accepted immutable result references.
#' @return NULL invisibly after required output acknowledgement.
#' @noRd
publishWorkflowPeptideResult <- function(owner, paths, result) {
    manager <- owner$manager
    entry <- workflowActionSpec(result$operation)
    record <- readDIAStateMetadata(manager, result$state_reference)
    envelope <- hydrateSavedWorkflowEntity(manager$.context$store, workflowPeptideProgressProduct(record$data))
    checkSavedWorkflowConsumerPayload("scientific_product", envelope)
    product <- envelope$value
    context_path <- file.path(paths$request_dir, "publication-context.rds")
    if (!file.exists(context_path)) {
        root <- file.path(paths$request_dir, "required-output")
        dir.create(root, showWarnings = FALSE)
        directories <- setupDirectories(root, entry$omic_type)
        directories$proteomics_dia_baseline_smoke <- directories$proteomics
        publication <- list(directories = directories,
            experiment_paths = directories$proteomics,
            omic_type = entry$omic_type, experiment_label = "dia_baseline_smoke")
        writeWorkflowExecutionMessage(context_path, publication)
    } else {
        publication <- readRDS(context_path)
        if (identical(sort(names(publication)),
            sort(c("experiment_label", "experiment_paths", "omic_type")))) {
            directories <- stats::setNames(list(publication$experiment_paths),
                paste(publication$omic_type, publication$experiment_label, sep = "_"))
            publication$directories <- directories
        } else if (is.list(publication) && is.list(publication$proteomics)) {
            publication <- list(directories = publication,
                experiment_paths = publication$proteomics,
                omic_type = entry$omic_type, experiment_label = "dia_baseline_smoke")
        }
    }
    if (!is.list(publication) || !is.list(publication$directories) ||
        !is.list(publication$experiment_paths) ||
        !identical(publication$omic_type, entry$omic_type) ||
        !is.character(publication$experiment_label) ||
        length(publication$experiment_label) != 1L || is.na(publication$experiment_label) ||
        !nzchar(publication$experiment_label) || hasMutableWorkflowEntityOwner(publication) ||
        as.numeric(object.size(publication)) > 65536) {
        cli::cli_abort("publication context is invalid", class = "workflow_execution_invalid")
    }
    directories <- publication$directories
    manager$.context$experiment_paths <- publication$experiment_paths
    # KIT-023: no global binding. The figure writer takes `project_dirs` as an argument and nothing in its chain
    # reads `config_list`; the controller is the application process, where a global would be shared by sessions.
    withr::local_options(multischolar.test_mode = FALSE)
    withr::local_envvar(MULTISCHOLAR_TEST_MODE = "false")
    grDevices::pdf(NULL)
    device <- grDevices::dev.cur()
    on.exit(if (device %in% grDevices::dev.list()) grDevices::dev.off(device))
    # The worker owns its filtering progress for this publication only (KIT-013): the product's
    # prior progress is seeded by the publisher, never read from or left in the global environment.
    arguments <- list(data = product$filteringSummary,
        step_name = workflowActionSpec(result$operation)$step_label,
        omic_type = publication$omic_type, experiment_label = publication$experiment_label,
        return_grid = TRUE, overwrite = TRUE, project_dirs = directories,
        progress_owner = new.env(parent = emptyenv()))
    workflowExecutionBoundary("publication")
    publishDIAQvalueOutput(manager, c(list(stateReference = result$state_reference), product),
        arguments, updateProteinFiltering)
    invisible(NULL)
}

#' Publish an accepted protein action's files and figures through the existing protein output path
#' (KIT-010): the required tables under the publication directories, then the filtering figures.
#' @param owner Leased controller.
#' @param paths Request-owned publication paths.
#' @param result Accepted immutable result references.
#' @return NULL invisibly after required output acknowledgement.
#' @noRd
publishWorkflowProteinResult <- function(owner, paths, result) {
    manager <- owner$manager
    entry <- workflowActionSpec(result$operation)
    reference <- result$state_reference
    product <- manager$.context$snapshot$output_products[[reference$entity_id]]
    if (is.null(product)) cli::cli_abort("protein output product is missing", class = "workflow_execution_invalid")
    envelope <- hydrateSavedWorkflowEntity(manager$.context$store, product)
    checkSavedWorkflowConsumerPayload("scientific_product", envelope)
    display <- envelope$value
    display$stateReference <- reference
    context_path <- file.path(paths$request_dir, "publication-context.rds")
    if (!file.exists(context_path)) {
        root <- file.path(paths$request_dir, "required-output")
        dir.create(root, showWarnings = FALSE)
        directories <- setupDirectories(root, entry$omic_type)
        directories$proteomics_dia_baseline_smoke <- directories$proteomics
        publication <- list(directories = directories, experiment_paths = directories$proteomics,
            omic_type = entry$omic_type, experiment_label = "dia_baseline_smoke")
        writeWorkflowExecutionMessage(context_path, publication)
    } else {
        publication <- readRDS(context_path)
        if (identical(sort(names(publication)), sort(c("experiment_label", "experiment_paths", "omic_type")))) {
            publication$directories <- stats::setNames(list(publication$experiment_paths),
                paste(publication$omic_type, publication$experiment_label, sep = "_"))
        } else if (is.list(publication) && is.list(publication$proteomics)) {
            publication <- list(directories = publication, experiment_paths = publication$proteomics,
                omic_type = entry$omic_type, experiment_label = "dia_baseline_smoke")
        }
    }
    if (!is.list(publication) || !is.list(publication$directories) || !is.list(publication$experiment_paths) ||
        !identical(publication$omic_type, entry$omic_type) || !is.character(publication$experiment_label) ||
        length(publication$experiment_label) != 1L || is.na(publication$experiment_label) ||
        !nzchar(publication$experiment_label) || hasMutableWorkflowEntityOwner(publication) ||
        as.numeric(object.size(publication)) > 65536) {
        cli::cli_abort("publication context is invalid", class = "workflow_execution_invalid")
    }
    manager$.context$experiment_paths <- publication$experiment_paths
    workflowExecutionBoundary("publication")
    publishDIAProteinTables(manager, display)
    # KIT-023: no global binding (see publishWorkflowPeptideResult).
    withr::local_options(multischolar.test_mode = FALSE)
    withr::local_envvar(MULTISCHOLAR_TEST_MODE = "false")
    grDevices::pdf(NULL)
    device <- grDevices::dev.cur()
    on.exit(if (device %in% grDevices::dev.list()) grDevices::dev.off(device))
    arguments <- list(data = display$filteringSummary, step_name = entry$step_label,
        omic_type = publication$omic_type, experiment_label = publication$experiment_label,
        return_grid = TRUE, overwrite = TRUE, project_dirs = publication$directories,
        progress_owner = new.env(parent = emptyenv()))
    publishDIAQvalueOutput(manager, display, arguments, updateProteinFiltering)
    invisible(NULL)
}

#' The diagnostic fixture's name: the one operation that is not a scientific declaration.
#' @noRd
workflowFixtureOperation <- "wfx_fixture"

#' Describe one registered operation (KIT-007: the registry entry, whose handler names never travel
#' in requests). An unregistered or malformed name refuses with `workflow_execution_unsupported`.
#' @param operation Internal operation identifier.
#' @return The registry entry: input role and class, handler names, limits, labels and variants.
#' @noRd
workflowActionSpec <- function(operation) workflowOperationEntry(operation)

#' Validate the registered scientific role without collecting assay data.
#' @param operation,manager,input_references Existing operation, owner and descriptors.
#' @return NULL invisibly; unsupported roles fail before worker dispatch.
#' @noRd
validateWorkflowActionInput <- function(operation, manager, input_references) {
    spec <- workflowActionSpec(operation)
    if (!is.list(input_references) || !identical(names(input_references), "state") ||
        !isSavedWorkflowReference(input_references$state))
        cli::cli_abort("unsupported input reference", class = "workflow_execution_invalid")
    record <- readDIAStateMetadata(manager, input_references$state)$data
    if (!isDIAComponentCheckpoint(record))
        cli::cli_abort("component checkpoint required", class = "workflow_execution_unsupported")
    if (spec$input == "peptide") {
        roles <- c("sample_id", "protein_id_column", "peptide_sequence_column")
        columns <- record$dimensions$columns
        valid <- identical(record$class, spec$input_class) &&
            all(c("peptide_data", "design_matrix") %in% names(record$components)) &&
            all(vapply(roles, function(role) {
                value <- record$slots[[role]]
                is.character(value) && length(value) == 1L && !is.na(value) &&
                    nzchar(value) && value %in% columns
            }, logical(1)))
        if (!valid) cli::cli_abort("peptide input roles are unavailable",
            class = "workflow_execution_unsupported")
    }
    invisible(NULL)
}

#' Freeze native-resolved arguments at the owned action boundary.
#' @param operation,manager,reference,parameters,policy Existing action inputs.
#' @return Bounded immutable arguments; no object or callback crosses the boundary.
#' @noRd
resolveWorkflowActionParameters <- function(operation, manager, reference, parameters, policy) {
    validateWorkflowActionParameters(operation, parameters, policy, resolved = FALSE)
    entry <- workflowActionSpec(operation)
    if (is.null(entry$resolve)) return(parameters)
    value <- get(entry$resolve, mode = "function")(manager, reference, parameters)
    validateWorkflowActionParameters(operation, value, policy)
    value
}

#' The sample filter's native resolver: the real method's parameter precedence on empty tables.
#' @noRd
resolveWorkflowSampleFilterParameters <- function(manager, reference, parameters) {
    prototype <- readDIAExecutionComponents(manager, reference, fields = character())
    # Invoke the real method's parameter resolver and call stack with empty tables.
    # Do not mirror its precedence or hydrate the scientific assay.
    arguments <- c(list(theObject = prototype, core_utilisation = NA), parameters)
    resolved <- do.call(filterMinNumPeptidesPerSample, arguments)
    native <- resolved@args$filterMinNumPeptidesPerSample
    list(peptides_per_sample_cutoff = native$peptides_per_sample_cutoff, inclusion_list = native$inclusion_list)
}

#' Complete the registered result only under the controller's existing acceptance lease.
#' @param owner,operation,value Leased controller, registered ID and immutable result.
#' @return Existing checkpoint/handoff result.
#' @noRd
completeWorkflowAction <- function(owner, operation, value, context = NULL, variant = NULL, values = NULL) {
    assertWorkflowController()
    if (is.null(owner$accepting)) cli::cli_abort("result acceptance lease required",
        class = "workflow_execution_busy")
    entry <- workflowActionSpec(operation)
    # KIT-025: the state this completion commits is the registered operation's, under the request's variant.
    setDIAStepProvenance(owner$manager, operation, variant %||% entry$default_variant %||% NA_character_,
        entry$parameters, entry$declaration$fixed, values = values)
    if (entry$completion == "fixture")
        return(commitDIAReferenceState(owner$manager, value$state_name, value$checkpoint,
            value$config, "guarded fixture computation"))
    workflow <- newWorkflowCompletionWorkflow(owner$manager, context)
    if (entry$completion == "protein") return(completeDIAProteinAction(workflow, value))
    if (entry$completion == "da") return(completeDIADaPublication(workflow, value))
    # KIT-014: the Q-value step keeps its own finishing step and product (`qvalue_progress`).
    if (identical(value$section, "srlQvalueProteotypicPeptideClean")) return(finishDIAReferenceQvalueStep(workflow, value))
    finishDIAReferencePeptideStep(workflow, value$checkpoint, value$state_name,
        value$description, value$config, value$section, value$counts,
        observed_inclusion_absence = if (isTRUE(entry$inclusion_absence))
            value$summary$inclusion_entries_absent_from_data else NULL)
}

#' Rank the repaired pilot on both widths without changing historical policies.
#' @param populations Primary/wide lists of native/physical/selection paired receipts.
#' @param expected_bindings Independently frozen per-population/per-arm bindings.
#' @return One minimax memory finalist, or a held/no-finalist decision; never activation.
#' @noRd
rankPerfRepairScreen <- function(populations, expected_bindings) {
    held <- list(status = "held", promotion_authorized = FALSE)
    if (!perfQualificationKeys(populations, c("primary", "wide")) ||
        !perfQualificationKeys(expected_bindings, c("primary", "wide"))) return(held)
    valid <- function(member, binding, rep) {
        is.list(member) && identical(member$status, "complete") &&
            identical(member$bindings, binding) && identical(member$rep, rep) &&
            all(vapply(c("science", "outputs", "measurement_valid", "integrity", "guard", "admission"),
                \(field) isTRUE(member[[field]]), logical(1))) &&
            is.numeric(member$peak_bytes) && length(member$peak_bytes) == 1L &&
            is.finite(member$peak_bytes) && member$peak_bytes > 0 && member$peak_bytes < 15.5 * 1024^3 &&
            is.numeric(member$elapsed_seconds) && length(member$elapsed_seconds) == 1L &&
            is.finite(member$elapsed_seconds) && member$elapsed_seconds > 0 && member$elapsed_seconds <= 300 &&
            identical(as.double(member$swap_bytes), 0) && identical(as.double(member$spill_bytes), 0)
    }
    medians <- list()
    eligible <- c("physical", "selection")
    requests <- character()
    for (population in names(populations)) {
        arms <- populations[[population]]
        bindings <- expected_bindings[[population]]
        if (!perfQualificationKeys(arms, c("native", "physical", "selection")) ||
            !perfQualificationKeys(bindings, names(arms)) ||
            any(!vapply(arms, \(arm) is.list(arm) && length(arm) == 2L, logical(1)))) return(held)
        for (arm in names(arms)) {
            binding <- bindings[[arm]]
            if (!is.character(binding) || length(binding) != 1L || is.na(binding) ||
                !grepl("^[a-f0-9]{64}$", binding)) return(held)
            for (member in arms[[arm]]) {
                id <- member$request_id
                if (!is.character(id) || length(id) != 1L || is.na(id) || !nzchar(id) || id %in% requests)
                    return(held)
                requests <- c(requests, id)
            }
            good <- all(vapply(seq_len(2L), \(i) valid(arms[[arm]][[i]], bindings[[arm]], i), logical(1)))
            if (!good && arm == "native") return(held)
            if (!good) eligible <- setdiff(eligible, arm)
            if (good) medians[[paste(population, arm)]] <- c(
                peak = stats::median(vapply(arms[[arm]], `[[`, numeric(1), "peak_bytes")),
                elapsed = stats::median(vapply(arms[[arm]], `[[`, numeric(1), "elapsed_seconds")))
        }
    }
    if (!length(eligible)) return(list(status = "no_finalist", promotion_authorized = FALSE))
    worst <- function(arm, metric) max(vapply(names(populations),
        \(population) medians[[paste(population, arm)]][[metric]], numeric(1)))
    peaks <- vapply(eligible, worst, numeric(1), metric = "peak")
    times <- vapply(eligible, worst, numeric(1), metric = "elapsed")
    tied <- which(peaks <= min(peaks) * 1.05)
    winner <- tied[order(times[tied], match(eligible[tied], c("physical", "selection")))][[1L]]
    list(status = "screen_finalist", finalist = eligible[[winner]],
        medians = medians, worst_peak = stats::setNames(peaks, eligible),
        worst_elapsed = stats::setNames(times, eligible), promotion_authorized = FALSE)
}

#' Frozen PERF policy, independent of historical schema-1 qualification.
#' @return Internal versioned thresholds and stage bounds.
#' @noRd
workflowPerfQualificationPolicy <- function() {
    list(schema = 2L, id = "perf.pilot.v1.2", peak_ratio = 0.80, elapsed_ratio = 1.20,
        screen_tie_ratio = 1.05, early_stop = 15.5 * 1024^3, memory_max = 16 * 1024^3,
        screen_seconds = 300, full_seconds = 1800, threads = 8L,
        engine_bytes = 12 * 1024^3, swap_bytes = 0, spill_bytes = 0)
}

#' Evaluate PERF screen feasibility separately from full paired benefit.
#' @param evidence Schema 2, stage, members/populations, immutable bindings and holds.
#' @param expected_bindings Independently frozen identities, including this policy's hash.
#' @param limits Optional tighter server-owned absolute memory/deadline caps.
#' @return Comparison decision and explicit non-promotion result.
#' @noRd
evaluatePerfWorkflowQualification <- function(evidence, expected_bindings, limits = NULL) {
    held <- list(status = "held", ratios = NULL, promotion_authorized = FALSE)
    policy <- workflowPerfQualificationPolicy()
    keys <- c("operation", "input", "parameters", "runtime", "boundary", "policy", "outputs",
        "execution_plan")
    bindings <- function(x) is.list(x) && length(x) == length(keys) &&
        !is.null(names(x)) && !anyNA(names(x)) && !anyDuplicated(names(x)) &&
        identical(sort(names(x)), sort(keys)) &&
        all(vapply(x, function(v) is.character(v) && length(v) == 1L && !is.na(v) &&
            grepl("^[a-f0-9]{64}$", v), logical(1)))
    flag <- function(x) is.logical(x) && length(x) == 1L && !is.na(x)
    positive <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0
    holds <- function(x) is.character(x) && !anyNA(x)
    if (!is.list(evidence) || !identical(evidence$schema, 2L) ||
        !identical(evidence$policy, policy$id) || !bindings(expected_bindings) ||
        !identical(expected_bindings$policy, digest::digest(policy, algo = "sha256")) ||
        !bindings(evidence$bindings) || !identical(evidence$bindings, expected_bindings) ||
        !flag(evidence$domain_supported) || !flag(evidence$runtime_eligible) ||
        !holds(evidence$open_holds) || !holds(evidence$release_holds) ||
        !is.character(evidence$goal) || length(evidence$goal) != 1L || is.na(evidence$goal) ||
        !evidence$goal %in% c("screen_feasibility", "benefit")) return(held)
    if (!identical(digest::digest(evidence$execution_plan, algo = "sha256"),
        expected_bindings$execution_plan) || !validatePerfQualificationMembers(evidence)) return(held)
    if (!evidence$domain_supported) return(utils::modifyList(held, list(status = "unsupported")))
    if (length(evidence$open_holds)) return(held)
    deadline <- if (evidence$goal == "screen_feasibility") policy$screen_seconds else policy$full_seconds
    memory_max <- policy$memory_max
    if (!is.null(limits)) {
        if (!is.list(limits) || !positive(limits$memory_max) || !positive(limits$deadline) ||
            limits$memory_max > memory_max || limits$deadline > deadline) return(held)
        memory_max <- limits$memory_max
        deadline <- limits$deadline
    }
    flags <- c("science", "outputs", "measurement_valid", "integrity", "guard", "admission")
    valid <- function(member) is.list(member) && identical(member$status, "complete") &&
        all(vapply(flags, function(name) isTRUE(member[[name]]), logical(1))) &&
        positive(member$peak_bytes) && positive(member$elapsed_seconds) &&
        member$peak_bytes < policy$early_stop && member$peak_bytes <= memory_max &&
        member$elapsed_seconds <= deadline && is.numeric(member$swap_bytes) &&
        identical(as.double(member$swap_bytes), 0) && is.numeric(member$spill_bytes) &&
        identical(as.double(member$spill_bytes), 0)
    if (evidence$goal == "screen_feasibility") {
        arms <- evidence$arms
        if (!perfQualificationKeys(arms, c("native", "physical", "selection")) ||
            !all(vapply(arms, function(x) is.list(x) && length(x) == 2L, logical(1)))) return(held)
        eligible <- c("physical", "selection")[vapply(arms[c("physical", "selection")],
            function(x) all(vapply(x, valid, logical(1))), logical(1))]
        if (!length(eligible)) return(utils::modifyList(held, list(status = "no_finalist")))
        peak <- vapply(arms[eligible], function(x) stats::median(vapply(x,
            function(m) m$peak_bytes, numeric(1))), numeric(1))
        elapsed <- vapply(arms[eligible], function(x) stats::median(vapply(x,
            function(m) m$elapsed_seconds, numeric(1))), numeric(1))
        tied <- which(peak <= min(peak) * policy$screen_tie_ratio)
        selected <- tied[order(elapsed[tied], match(eligible[tied], c("physical", "selection")))][[1L]]
        return(list(status = "screen_finalist", finalist = eligible[[selected]],
            medians = list(peak = peak, elapsed = elapsed), ratios = NULL,
            promotion_authorized = FALSE))
    }
    populations <- evidence$populations
    if (!perfQualificationKeys(populations, c("primary", "wide")) ||
        !all(vapply(populations, function(x) is.list(x) && length(x) == 2L, logical(1)))) return(held)
    ratios <- list()
    for (population in c("primary", "wide")) {
        ratios[[population]] <- list()
        for (i in seq_len(2L)) {
            pair <- populations[[population]][[i]]
            if (!is.list(pair) || !valid(pair$control)) return(held)
            if (!valid(pair$candidate)) return(utils::modifyList(held, list(status = "rejected")))
            ratios[[population]][[i]] <- list(peak = pair$candidate$peak_bytes / pair$control$peak_bytes,
                elapsed = pair$candidate$elapsed_seconds / pair$control$elapsed_seconds)
        }
    }
    passed <- unlist(lapply(ratios, function(pairs) vapply(pairs, function(pair)
        pair$peak <= policy$peak_ratio && pair$elapsed <= policy$elapsed_ratio, logical(1))))
    comparison <- if (all(passed)) "benefit" else if (any(passed)) "inconclusive" else "rejected"
    list(status = if (!evidence$runtime_eligible || length(evidence$release_holds)) "held" else comparison,
        comparison_status = comparison, ratios = ratios, promotion_authorized = FALSE)
}

#' Require exact named collection keys without dropping missing names.
#' @param value,keys Collection and its fixed key set.
#' @return TRUE only for the declared collection shape.
#' @noRd
perfQualificationKeys <- function(value, keys) {
    is.list(value) && length(value) == length(keys) && !is.null(names(value)) &&
        !anyNA(names(value)) && !anyDuplicated(names(value)) &&
        all(nzchar(names(value))) && setequal(names(value), keys)
}

#' Bind every PERF member to the independently frozen execution plan.
#' @param evidence Schema-2 evidence; its plan hash was independently checked.
#' @return FALSE on missing, mixed, duplicated or mismatched member identities.
#' @noRd
validatePerfQualificationMembers <- function(evidence) {
    plan <- evidence$execution_plan
    if (!perfQualificationKeys(plan, c("schema", "goal", "finalist", "screen_decision", "members")) ||
        !identical(plan$schema, 1L) || !identical(plan$goal, evidence$goal)) return(FALSE)
    scalar <- function(x) is.character(x) && length(x) == 1L && !is.na(x) &&
        nzchar(x) && nchar(x, type = "bytes") <= 128L
    hash <- function(x) scalar(x) && grepl("^[a-f0-9]{64}$", x)
    pair_list <- function(x) is.list(x) && length(x) == 2L && is.null(names(x))
    observed <- expected_arm <- expected_population <- expected_rep <- list()
    if (evidence$goal == "screen_feasibility") {
        if (!is.null(plan$finalist) || !is.null(plan$screen_decision) ||
            !perfQualificationKeys(evidence$arms, c("native", "physical", "selection"))) return(FALSE)
        for (arm in c("native", "physical", "selection")) {
            if (!pair_list(evidence$arms[[arm]])) return(FALSE)
            for (rep in seq_len(2L)) {
                slot <- paste("screen", arm, rep, sep = ".")
                observed[[slot]] <- evidence$arms[[arm]][[rep]]
                expected_arm[[slot]] <- arm
                expected_population[[slot]] <- "screen"
                expected_rep[[slot]] <- rep
            }
        }
    } else {
        if (!scalar(plan$finalist) || !plan$finalist %in% c("physical", "selection") ||
            !hash(plan$screen_decision) ||
            !perfQualificationKeys(evidence$populations, c("primary", "wide"))) return(FALSE)
        for (population in c("primary", "wide")) {
            if (!pair_list(evidence$populations[[population]])) return(FALSE)
            for (rep in seq_len(2L)) {
                pair <- evidence$populations[[population]][[rep]]
                if (!perfQualificationKeys(pair, c("control", "candidate"))) return(FALSE)
                for (arm in c("control", "candidate")) {
                    slot <- paste(population, arm, rep, sep = ".")
                    observed[[slot]] <- pair[[arm]]
                    expected_arm[[slot]] <- if (arm == "control") "native" else plan$finalist
                    expected_population[[slot]] <- population
                    expected_rep[[slot]] <- rep
                }
            }
        }
    }
    if (!perfQualificationKeys(plan$members, names(expected_arm)) ||
        !perfQualificationKeys(observed, names(expected_arm))) return(FALSE)
    identity_keys <- c("run_id", "population", "arm", "repetition", "input", "parameters",
        "runtime", "build", "boundary", "output_contract")
    hashes <- c("input", "parameters", "runtime", "build", "boundary", "output_contract")
    for (slot in names(expected_arm)) {
        wanted <- plan$members[[slot]]
        member <- observed[[slot]]
        if (!perfQualificationKeys(wanted, identity_keys) || !scalar(wanted$run_id) ||
            !all(vapply(wanted[hashes], hash, logical(1))) ||
            !identical(wanted$arm, expected_arm[[slot]]) ||
            !identical(wanted$population, expected_population[[slot]]) ||
            !identical(wanted$repetition, expected_rep[[slot]]) || !is.list(member) ||
            is.null(names(member)) || anyNA(names(member)) || anyDuplicated(names(member)) ||
            !all(c(identity_keys, "receipt_id") %in% names(member)) || !hash(member$receipt_id) ||
            !identical(member[identity_keys], wanted[identity_keys])) return(FALSE)
    }
    values <- function(field) vapply(observed, function(x) x[[field]], character(1))
    if (anyDuplicated(values("run_id")) || anyDuplicated(values("receipt_id"))) return(FALSE)
    # Same scientific inputs/parameters/output contract for every arm and repeat
    # of one population; the frozen build/runtime stays fixed for each arm.
    for (population in unique(values("population"))) for (field in c("input", "parameters", "boundary", "output_contract")) {
        if (length(unique(values(field)[values("population") == population])) != 1L) return(FALSE)
    }
    for (arm in unique(values("arm"))) for (field in c("build", "runtime")) {
        if (length(unique(values(field)[values("arm") == arm])) != 1L) return(FALSE)
    }
    if (evidence$goal == "benefit" && length(unique(values("input"))) != 2L) return(FALSE)
    TRUE
}
