# Internal WFX execution boundary. Scientific adapters are deliberately finite.

#' Evaluate frozen paired qualification evidence without executing an operation.
#' @param evidence Bounded typed qualification record (schema 1).
#' @param expected_bindings Independently frozen operation/input/runtime/boundary identities.
#' @param limits Server-owned absolute caps; NULL uses the common qualification maximum.
#' @return One deterministic decision; malformed or incomplete evidence is held.
#' @noRd
evaluateWorkflowQualification <- function(evidence, expected_bindings, limits = NULL) {
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
        operation = pending$operation %||% "wfx_fixture")
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
        if (!identical(current$entity_id, proposed$entity_id))
            cli::cli_abort("head recovery conflicts with a later mutation", class = "workflow_execution_stale")
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
    value <- suppressWarnings(as.numeric(strsplit(trimws(line), "[[:space:]]+")[[1L]][[2L]])) * 1024
    if (length(value) != 1L || !is.finite(value) || value < 0) cli::cli_abort("host memory unavailable")
    value
}

#' Return the server-owned policy; small diagnostics never qualify scientific inputs.
#' @param diagnostic Explicit internal lifecycle-fixture mode.
#' @return Fixed resource envelope, never populated from GUI input.
#' @noRd
workflowExecutionPolicy <- function(diagnostic = FALSE) {
    if (diagnostic) return(list(schema = 1L, diagnostic = TRUE, threads = 8L,
        memory_max = 3 * 1024^3, early_stop = 2.9 * 1024^3, swap_max = 0,
        growth_reserve = 128 * 1024^2, startup_headroom = 2 * 1024^3,
        host_floor = 2 * 1024^3, engine_limit = 12 * 1024^3, deadline = 90L,
        monitor_interval = 0.05))
    list(schema = 1L, diagnostic = FALSE, threads = 8L, memory_max = 16 * 1024^3,
        early_stop = 15.5 * 1024^3, swap_max = 0, growth_reserve = 1024^3,
        startup_headroom = 4 * 1024^3, host_floor = 2 * 1024^3,
        engine_limit = 12 * 1024^3, deadline = 1800L, monitor_interval = 0.05)
}

#' Managed application resource owner; no scientific state is stored here.
#' @noRd
WorkflowExecutionApplication <- R6::R6Class("WorkflowExecutionApplication", cloneable = FALSE,
    public = list(config = NULL, cgroup = NULL, controllers = list(), active = NULL,
        sessions = 0L, pid = NULL, closed = FALSE))

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
    if (!isTRUE(actual$available) || !identical(as.numeric(actual$MainPID), as.numeric(Sys.getpid())) ||
        !startsWith(actual$ControlGroup, paste0(config$control_group, "/")) ||
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
    con <- file(file.path(application$cgroup, "memory.peak"), open = "r+")
    ready <- FALSE
    on.exit(if (!ready) close(con))
    writeLines("0", con)
    flush(con)
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
launchWorkflowApplication <- function(root, diagnostic = FALSE, mode = "gui", test_root = NULL, port = 3838L) {
    if (!identical(Sys.info()[["sysname"]], "Linux") || !dir.exists(root) ||
        !mode %in% c("gui", "protocol", "resources", "adapters") || mode %in% c("protocol", "resources") && !isTRUE(diagnostic))
        cli::cli_abort("unsupported managed application launch", class = "workflow_execution_unsupported")
    policy <- workflowExecutionPolicy(diagnostic)
    if (workflowHostAvailable() < policy$startup_headroom)
        cli::cli_abort("insufficient host headroom", class = "workflow_execution_refused")
    workflowSystemd(c("show-environment"))
    id <- substr(digest::digest(list(Sys.time(), Sys.getpid(), tempfile()), algo = "sha256"), 1L, 16L)
    stem <- paste0("multischolarwfx", id)
    slice <- paste0(stem, ".slice")
    config <- list(schema = 1L, root = normalizePath(root), stem = stem, slice = slice,
        controller_unit = paste0(stem, "-controller.service"), policy = policy,
        build = workflowExecutionBuildIdentity(), library = dirname(system.file(package = "MultiScholaR")),
        mode = mode, port = as.integer(port), test_root = test_root)
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
        config$test_hashes <- stats::setNames(vapply(inputs, function(path)
            digest::digest(file = path, algo = "sha256"), character(1)), inputs)
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
    workflowLaunchService(config, config$controller_unit, "application", config_path,
        digest::digest(file = config_path, algo = "sha256"), deadline = if (mode == "adapters") { if (diagnostic) 300L else 1800L } else
            if (diagnostic) 600L else 86400L)
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
    on.exit({
        for (owner in app$controllers) closeWorkflowExecutionController(owner)
        app$closed <- TRUE
    })
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
    shiny::runApp(run_app(workflow_execution = app), host = "127.0.0.1",
        port = config$port, launch.browser = FALSE)
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
            workflowSystemd(c("stop", request$worker_unit), required = FALSE)
            writeWorkflowExecutionMessage(file.path(request$directory, "guard-stop.rds"),
                list(schema = 1L, reason = stop_reason, current = current, host_available = host))
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
    if (!isTRUE(request$config$policy$diagnostic) || request$operation != "wfx_fixture")
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
    list(checkpoint = record$data, config = record$config, state_name = "wfx_fixture",
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
        value <- switch(request$operation, wfx_fixture = computeWorkflowFixture(request),
            precursor_rollup = computeWorkflowPeptideAction(request),
            minimum_peptides_per_sample = computeWorkflowPeptideAction(request),
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
#' @return Request ID immediately after monitor launch; pollWorkflowAction advances it.
#' @noRd
runWorkflowAction <- function(operation, input_references, parameters, resource_policy) {
    entered <- proc.time()[["elapsed"]]
    owner <- resource_policy
    if (!inherits(owner, "WorkflowExecutionController") || owner$closed) cli::cli_abort("controller required")
    app <- owner$application
    peak <- openWorkflowPeakWindow(app)
    transferred <- FALSE
    on.exit(if (!transferred) close(peak), add = TRUE)
    checkWorkflowApplication(app)
    assertDIAExportReady(owner$manager)
    if (!is.null(app$active) || !is.null(owner$job)) cli::cli_abort("another heavy action is active", class = "workflow_execution_busy")
    validateWorkflowActionParameters(operation, parameters, app$config$policy)
    if (!is.list(input_references) || !identical(names(input_references), "state") ||
        !isSavedWorkflowReference(input_references$state))
        cli::cli_abort("unsupported input reference", class = "workflow_execution_invalid")
    if (workflowHostAvailable() < app$config$policy$startup_headroom)
        cli::cli_abort("insufficient host headroom", class = "workflow_execution_refused")
    admission <- filelock::lock(file.path(app$config$root, ".heavy-action.lock"), timeout = 0)
    if (is.null(admission)) cli::cli_abort("application already has a heavy action", class = "workflow_execution_busy")
    ready <- FALSE
    on.exit(if (!ready) {
        filelock::unlock(admission)
        owner$job <- NULL
        app$active <- NULL
    })
    job <- new.env(parent = emptyenv())
    job$started <- entered
    job$peak <- peak
    job$accounting <- workflowCgroupAccounting(app$cgroup)
    store <- owner$manager$.context$store
    validateSavedWorkflowReference(store, input_references$state, payload = FALSE)
    if (!identical(input_references$state, getDIAStateReference(owner$manager)))
        cli::cli_abort("action input is not the selected immutable checkpoint", class = "workflow_execution_stale")
    readDIAStateMetadata(owner$manager, input_references$state)
    id <- paste0("wfx-", substr(digest::digest(list(Sys.time(), tempfile(), owner$control$generation), algo = "sha256"), 1, 32))
    paths <- workflowExportPaths(store$directory, id)
    token <- list(request_id = id, incarnation = store$identity$incarnation,
        input_ids = vapply(input_references, `[[`, character(1), "entity_id"),
        expected_head = workflowExecutionHead(store)$entity_id, generation = owner$control$generation,
        controller_id = owner$control$controller_id, parameters = digest::digest(parameters, algo = "sha256"),
        execution = workflowExecutionBuildIdentity())
    request <- list(schema = 1L, build = token$execution, operation = operation, input_references = input_references,
        parameters = parameters, token = token, owner = store$identity, project_root = dirname(store$directory),
        directory = paths$request_dir, config = app$config)
    job$request <- request
    job$paths <- paths
    job$admission <- admission
    job$stage <- "admitting"
    job$worker_unit <- paste0(app$config$stem, "-w", substr(token$request_id, 5, 16), ".service")
    job$monitor_unit <- paste0(app$config$stem, "-m", substr(token$request_id, 5, 16), ".service")
    job$worker_request <- file.path(paths$request_dir, "request.rds")
    writeWorkflowExecutionMessage(job$worker_request, request)
    monitor <- list(schema = 1L, build = request$build, config = app$config,
        directory = paths$request_dir, worker_unit = job$worker_unit)
    monitor_path <- file.path(paths$request_dir, "monitor.rds")
    writeWorkflowExecutionMessage(monitor_path, monitor)
    writeWorkflowExportJournal(paths, list(schema = "multischolar.workflow_action.v1", request_id = id,
        status = "admitting", expected_head = token$expected_head, token_hash = digest::digest(token, algo = "sha256")))
    owner$control$active_request <- id
    writeWorkflowExecutionControl(store, owner$control)
    owner$job <- job
    app$active <- id
    # Transfer ownership before launch: a failed acknowledgement may still have created M.
    ready <- TRUE
    transferred <- TRUE
    launched <- tryCatch({
        workflowLaunchService(app$config, job$monitor_unit, "monitor", monitor_path,
            digest::digest(file = monitor_path, algo = "sha256"), binds = app$config$controller_unit,
            deadline = app$config$policy$deadline + 60L)
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
    guard <- file.path(job$paths$request_dir, "guard-stop.rds")
    if (file.exists(guard)) return(finishWorkflowExecutionJob(owner, readRDS(guard)$reason))
    monitor <- workflowServiceProperties(job$monitor_unit)
    if (!isTRUE(monitor$available) || !monitor$SubState %in% c("running", "start"))
        return(finishWorkflowExecutionJob(owner, "monitor_lost"))
    if (job$stage == "admitting") {
        if (!file.exists(file.path(job$paths$request_dir, "monitor-ready.rds"))) {
            if (proc.time()[["elapsed"]] - job$started > 30) return(finishWorkflowExecutionJob(owner, "monitor_start_timeout"))
            return(list(status = "admitting"))
        }
        current <- workflowMemoryCounter(file.path(app$cgroup, "memory.current"))
        reservation <- current + app$config$policy$growth_reserve
        cap <- app$config$policy$memory_max - reservation
        # This finite fixture estimate includes R/package startup and bounded reader overhead.
        if (cap < 768 * 1024^2) return(finishWorkflowExecutionJob(owner, "admission_refused"))
        writeWorkflowExecutionMessage(file.path(job$paths$request_dir, "admitted.rds"),
            list(schema = 1L, reservation = reservation, worker_cap = cap))
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
    result_path <- file.path(job$paths$request_dir, "worker-result.rds")
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
    if (identical(job$request$operation, "wfx_fixture")) {
        commitDIAReferenceState(owner$manager, value$state_name, value$checkpoint, value$config,
            "guarded fixture computation")
    } else {
        workflow <- new.env(parent = emptyenv())
        workflow$state_manager <- owner$manager
        workflow$qc_params <- list(peptide_qc = list())
        finishDIAReferencePeptideStep(workflow, value$checkpoint, value$state_name,
            value$description, value$config, value$section, value$counts,
            observed_inclusion_absence = if (job$request$operation == "minimum_peptides_per_sample")
                value$summary$inclusion_entries_absent_from_data else NULL)
    }
    owner$accepting <- NULL
    journal <- readWorkflowExecutionJournal(job$paths)
    journal$status <- "publication_pending"
    writeWorkflowExportJournal(job$paths, journal)
    published <- tryCatch({ publishWorkflowFixtureResult(owner, job$paths); TRUE },
        error = function(error) FALSE)
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
    if (!identical(result$operation %||% "wfx_fixture", "wfx_fixture")) {
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
    if (!valid_store) workflowSystemd(c("stop", job$monitor_unit), required = FALSE)
    worker <- workflowServiceProperties(job$worker_unit)
    workflowSystemd(c("stop", job$worker_unit), required = FALSE)
    if (valid_store && dir.exists(job$paths$request_dir))
        writeWorkflowExecutionMessage(file.path(job$paths$request_dir, "monitor-done.rds"), list(schema = 1L))
    deadline <- proc.time()[["elapsed"]] + 2
    monitor_path <- file.path(job$paths$request_dir, "monitor-result.rds")
    while (valid_store && !file.exists(monitor_path) && proc.time()[["elapsed"]] < deadline) Sys.sleep(0.02)
    workflowSystemd(c("stop", job$monitor_unit), required = FALSE)
    worker_empty <- workflowServiceEmpty(owner$application, job$worker_unit)
    monitor_empty <- workflowServiceEmpty(owner$application, job$monitor_unit)
    if (!worker_empty || !monitor_empty) {
        job$stage <- "cleanup_failed"
        job$terminal_status <- status
        return(list(status = "cleanup_failed", request_id = job$request$token$request_id))
    }
    on.exit({
        filelock::unlock(job$admission)
        owner$application$active <- NULL
        owner$job <- NULL
    })
    if (valid_store) for (role in c("worker", "monitor")) {
        scratch <- file.path(job$paths$request_dir, paste0(role, "-scratch"))
        if (dir.exists(scratch) && !nzchar(Sys.readlink(scratch))) unlink(scratch, recursive = TRUE)
    }
    peak <- closeWorkflowPeakWindow(job$peak)
    accounting <- workflowCgroupAccounting(owner$application$cgroup)
    delta <- Map(`-`, accounting, job$accounting)
    files <- if (valid_store) list.files(job$paths$request_dir, recursive = TRUE, full.names = TRUE) else character()
    staged_bytes <- if (valid_store) sum(file.info(files)$size) else NA_real_
    monitor <- if (valid_store && file.exists(monitor_path)) readRDS(monitor_path) else NULL
    receipt <- list(schema = 1L, request_id = job$request$token$request_id, status = status,
        token = job$request$token, error = job$error, elapsed_seconds = proc.time()[["elapsed"]] - job$started,
        peak_bytes = peak, whole_service_peak_bytes = workflowMemoryCounter(file.path(owner$application$cgroup, "memory.peak")),
        monitor = monitor, worker = worker, threads = 8L, browser_memory_included = FALSE,
        accounting_delta = delta, staged_bytes_at_teardown = staged_bytes,
        spill_attempts = if (job$request$operation == "wfx_fixture") 0L else NA_integer_,
        spill_basis = if (job$request$operation == "wfx_fixture")
            "identity fixture uses no database or spill-capable algorithm" else
            "engine spill disabled; attempted spill is not independently counted",
        measurement_valid = !is.null(monitor) && monitor$max_observation_interval <= 0.1,
        cleanup = list(worker_stopped = worker_empty, monitor_stopped = monitor_empty))
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
workflowExecutionProcess <- function() {
    if (!identical(Sys.info()[["sysname"]], "Linux")) return(list(pid = Sys.getpid(), platform = Sys.info()[["sysname"]]))
    stat <- strsplit(sub("^.*[)] ", "", readLines("/proc/self/stat", warn = FALSE)), " ", fixed = TRUE)[[1L]]
    list(pid = Sys.getpid(), start = stat[[20L]], boot = readLines("/proc/sys/kernel/random/boot_id", warn = FALSE))
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
    runWorkflowAction("wfx_fixture", list(state = getDIAStateReference(manager)), list(mode = "identity"), owner)
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
#' @return NULL invisibly or explicit refusal before dispatch.
#' @noRd
validateWorkflowActionParameters <- function(operation, parameters, policy) {
    if (!is.character(operation) || length(operation) != 1L || is.na(operation) ||
        !operation %in% c("wfx_fixture", "precursor_rollup", "minimum_peptides_per_sample"))
        cli::cli_abort("operation adapter is not qualified", class = "workflow_execution_unsupported")
    # The frozen 4,000-sample filter includes duplicate/whitespace-bearing inclusion
    # arguments. Preserve those native argument bytes under a finite operation-specific cap.
    parameter_limit <- if (operation == "minimum_peptides_per_sample") 262144 else 65536
    if (!is.list(parameters) || hasMutableWorkflowEntityOwner(parameters) ||
        as.numeric(object.size(parameters)) > parameter_limit || anyDuplicated(names(parameters)))
        cli::cli_abort("unsupported execution parameters", class = "workflow_execution_invalid")
    if (operation == "precursor_rollup") {
        if (length(parameters)) cli::cli_abort("P14 has no caller-selectable execution settings")
        diaRollupExtensionInfo()
    } else if (operation == "minimum_peptides_per_sample") {
        if (is.null(names(parameters)) || !setequal(names(parameters),
            c("peptides_per_sample_cutoff", "inclusion_list"))) cli::cli_abort("sample filter parameters required")
        k <- parameters$peptides_per_sample_cutoff
        if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k < 1 ||
            k > .Machine$integer.max || k != floor(k)) cli::cli_abort("cutoff must be a positive integer")
        if (!is.null(parameters$inclusion_list) && !is.character(parameters$inclusion_list))
            cli::cli_abort("inclusion entries must be character")
    } else {
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
    }
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
    value <- switch(request$operation,
        precursor_rollup = computeDIAReferencePrecursorRollup(manager, request$input_references$state,
            worker_cap = request$worker_cap),
        minimum_peptides_per_sample = computeDIAReferenceSampleFilter(manager,
            request$input_references$state, request$parameters, worker_cap = request$worker_cap),
        cli::cli_abort("scientific adapter unavailable"))
    value$operation <- request$operation
    value
}

#' Publish accepted peptide QC figures using existing checkpoint/output recovery.
#' @param owner Leased controller.
#' @param paths Request-owned publication paths.
#' @param result Accepted immutable result references.
#' @return NULL invisibly after required output acknowledgement.
#' @noRd
publishWorkflowPeptideResult <- function(owner, paths, result) {
    manager <- owner$manager
    record <- readDIAStateMetadata(manager, result$state_reference)
    envelope <- hydrateSavedWorkflowEntity(manager$.context$store, record$data$products$qc_progress)
    checkSavedWorkflowConsumerPayload("scientific_product", envelope)
    product <- envelope$value
    context_path <- file.path(paths$request_dir, "publication-context.rds")
    if (!file.exists(context_path)) {
        root <- file.path(paths$request_dir, "required-output")
        dir.create(root, showWarnings = FALSE)
        directories <- setupDirectories(root, "proteomics")
        directories$proteomics_dia_baseline_smoke <- directories$proteomics
        writeWorkflowExecutionMessage(context_path, directories)
    }
    directories <- readRDS(context_path)
    manager$.context$experiment_paths <- directories$proteomics
    rlang::local_bindings(config_list = record$config, project_dirs = directories,
        filtering_progress = methods::new("FilteringProgress"), .env = .GlobalEnv)
    withr::local_options(multischolar.test_mode = FALSE)
    withr::local_envvar(MULTISCHOLAR_TEST_MODE = "false")
    grDevices::pdf(NULL)
    device <- grDevices::dev.cur()
    on.exit(if (device %in% grDevices::dev.list()) grDevices::dev.off(device))
    arguments <- list(data = product$filteringSummary,
        step_name = if (result$operation == "precursor_rollup") "4_precursor_rollup" else "6_sample_filtered",
        omic_type = "proteomics", experiment_label = "dia_baseline_smoke",
        return_grid = TRUE, overwrite = TRUE, project_dirs = directories)
    workflowExecutionBoundary("publication")
    publishDIAQvalueOutput(manager, c(list(stateReference = result$state_reference), product),
        arguments, updateProteinFiltering)
    invisible(NULL)
}
