# KIT stage coordinator (modelled on workflow_perf106_repair.R; CORE-pinned files untouched).
# Usage: Rscript workflow_kit_coordinator.R <stage> <library> <root> [--validate-only]
# The coordinator is a metadata-only launcher outside the managed slice: it launches the
# application slice, samples its memory, drives the browser and verifies afterwards.
args <- commandArgs(trailingOnly = TRUE)
validate_only <- length(args) == 4L && identical(args[[4L]], "--validate-only")
if (validate_only) args <- args[1:3]
stopifnot(length(args) == 3L, grepl("^kit-(lifecycle-(complete|fail)|journey)$", args[[1L]]))
stage <- args[[1L]]
library_path <- normalizePath(args[[2L]])
root <- args[[3L]]
repo <- normalizePath(".")
.libPaths(c(library_path, .libPaths()))
loadNamespace("MultiScholaR")
api <- function(name) get(name, asNamespace("MultiScholaR"))
dir.create(root, recursive = TRUE, showWarnings = FALSE)
stage_root <- if (validate_only) file.path(tempfile("kit-validate-"), stage) else file.path(root, stage)
if (validate_only) dir.create(dirname(stage_root))
stopifnot(!file.exists(stage_root), dir.create(stage_root))
started <- Sys.time()
receipt_path <- file.path(stage_root, "stage.rds")
receipt <- list(stage = stage, status = "starting", started = started, cases = list())
saveRDS(receipt, receipt_path)
snapshot <- file.path(stage_root, "source")
files <- c("tools/profiling/workflow_kit.R", "tools/profiling/workflow_kit_driver.R",
    "tools/profiling/workflow_kit_coordinator.R", "tools/profiling/run_workflow_kit.sh",
    "tools/profiling/workflow_checkpoint_dia_driver.R",
    "docs/qa/workflow-execution/kit/bindings.json", "docs/qa/workflow-execution/kit/protocol.md",
    "docs/qa/workflow-execution/kit/kit-001/step-map.json")
for (file in files) {
    target <- file.path(snapshot, file)
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    stopifnot(file.copy(file.path(repo, file), target))
}
receipt$source_hashes <- stats::setNames(vapply(file.path(snapshot, files),
    \(file) digest::digest(file = file, algo = "sha256"), character(1)), files)
helpers <- new.env(parent = globalenv())
sys.source(file.path(snapshot, "tools/profiling/workflow_checkpoint_dia_driver.R"), helpers)
sys.source(file.path(snapshot, "tools/profiling/workflow_kit.R"), helpers)
sys.source(file.path(snapshot, "tools/profiling/workflow_kit_driver.R"), helpers)
bindings_path <- file.path(snapshot, "docs/qa/workflow-execution/kit/bindings.json")
bindings <- helpers$kitBindings(bindings_path)
cases <- helpers$kitCases(bindings, stage)
member_seconds <- helpers$kitMemberSeconds(stage, bindings)
receipt$input_bindings_sha256 <- digest::digest(file = bindings_path, algo = "sha256")
receipt$build <- api("workflowExecutionBuildIdentity")()
hash <- function(path) digest::digest(file = path, algo = "sha256")
stopifnot(identical(normalizePath(root), normalizePath(bindings$root)),
    identical(hash(file.path(snapshot, "docs/qa/workflow-execution/kit/protocol.md")), bindings$protocol_sha256),
    identical(hash(file.path(snapshot, "docs/qa/workflow-execution/kit/kit-001/step-map.json")), bindings$step_map_sha256),
    identical(unname(unlist(bindings$files)), vapply(file.path(repo, names(bindings$files)), hash, character(1), USE.NAMES = FALSE)),
    normalizePath(library_path) %in% normalizePath(vapply(cases, `[[`, character(1), "library")))
for (case in cases) {
    stopifnot(dir.exists(case$library))
    corpus <- bindings$corpora[[case$corpus]]
    for (name in names(corpus$files)) stopifnot(identical(hash(corpus$files[[name]]$path), corpus$files[[name]]$sha256))
}
hostReady <- function() {
    units <- system2("systemctl", c("--user", "list-units", "--state=active", "--no-legend", "--plain",
        shQuote("multischolar*")), stdout = TRUE)
    profile_path <- "/sys/firmware/acpi/platform_profile"
    cpu <- function() as.numeric(strsplit(readLines("/proc/stat", n = 1L), " +")[[1L]][-1L])
    before <- cpu()
    Sys.sleep(2)
    ticks <- cpu() - before
    list(available = api("workflowHostAvailable")(), units = units[nzchar(trimws(units))],
        profile = if (file.exists(profile_path)) readLines(profile_path, n = 1L) else NA_character_,
        cpu_busy = 1 - sum(ticks[4:5]) / sum(ticks[1:8]))
}
receipt$host_at_start <- hostReady()
if (validate_only) {
    refusal <- helpers$kitSetupRefusal(receipt$host_at_start)
    if (!is.null(refusal)) {
        cat("SETUP_REFUSED", stage, refusal, "\n")
        unlink(dirname(stage_root), recursive = TRUE)
        quit(save = "no", status = 76L)
    }
    cat("VALIDATED", stage, length(cases), "case(s); host", format(receipt$host_at_start$available / 1024^3,
        digits = 3), "GiB available;", length(receipt$host_at_start$units), "active managed unit(s)\n")
    unlink(dirname(stage_root), recursive = TRUE)
    quit(save = "no", status = 0L)
}
clock <- helpers$kitClock
sliceState <- function(slice) {
    state <- system2("systemctl", c("--user", "show", shQuote(slice), "--property=ActiveState", "--value"), stdout = TRUE)
    if (!length(state)) "unknown" else trimws(state[[1L]])
}
active_slice <- NULL
sampler <- NULL
tryCatch({
    for (case in cases) {
        host <- hostReady()
        refusal <- helpers$kitSetupRefusal(host)
        if (!is.null(refusal)) {
            receipt$status <- "setup_refused"
            receipt$refusal <- list(member = case$id, reason = refusal, host = host)
            saveRDS(receipt, receipt_path)
            stop("KIT setup refused: ", refusal)
        }
        owned <- file.path(stage_root, case$id)
        stopifnot(dir.create(owned))
        application_root <- file.path(owned, "application")
        dir.create(application_root)
        member <- bindings$journey$member
        stopifnot(identical(member$programme, "kit-journey-v1"))
        saveRDS(list(schema = 1L, programme = member$programme, deadline = as.numeric(member$deadline),
            guard_stop_test = case$expect_failure), file.path(application_root, "gui-member.rds"),
            compress = FALSE, version = 3)
        saveRDS(c(list(schema = 1L, stage = stage, deadline = member_seconds), case), file.path(owned, "timed-member.rds"))
        port <- httpuv::randomPort()
        # Launch from a fresh process on the case library; the coordinator never swaps namespaces.
        expression <- paste(
            'a <- commandArgs(TRUE); .libPaths(c(a[[1L]], .libPaths()));',
            'ns <- loadNamespace("MultiScholaR");',
            'stamp <- as.numeric(strsplit(readLines("/proc/uptime", n = 1L), " ")[[1L]][[1L]]);',
            'saveRDS(stamp, file.path(a[[2L]], "launcher-start.rds"));',
            'get("launchWorkflowApplication", ns)(a[[3L]], diagnostic = FALSE, mode = "gui", port = as.integer(a[[4L]]))')
        status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", "-e", shQuote(expression),
            shQuote(c(case$library, owned, application_root, port))))
        stopifnot(identical(status, 0L))
        config <- readRDS(file.path(application_root, "application.rds"))
        operation_started <- readRDS(file.path(owned, "launcher-start.rds"))
        deadline_uptime <- operation_started + member_seconds
        stopifnot(identical(config$build, case$expected_build), identical(config$mode, "gui"),
            isTRUE(config$journey$native_sessions))
        cat("START", case$id, config$controller_unit, "\n")
        active_slice <- config$slice
        group <- file.path("/sys/fs/cgroup", sub("^/", "", config$control_group))
        trace_path <- file.path(owned, "memory-trace.tsv")
        stop_file <- file.path(owned, "stop-sampling")
        sampler_log <- file.path(owned, "sampler.log")
        system2("bash", shQuote(c(file.path(snapshot, "tools/profiling/run_workflow_kit.sh"), "--sample", group,
            trace_path, sprintf("%.3f", deadline_uptime), config$slice, stop_file, file.path(owned, "guard-stop.txt"))),
            stdout = sampler_log, stderr = sampler_log, wait = FALSE)
        sampler <- stop_file
        corpus_spec <- bindings$corpora[[case$corpus]]
        corpus <- list(report = corpus_spec$files[["report.tsv"]]$path, fasta = corpus_spec$files[["proteins.fasta"]]$path,
            design = corpus_spec$files[["design.tsv"]]$path)
        journey <- bindings$journey
        journey$corpus_key <- case$corpus
        fault <- if (case$expect_failure) bindings$lifecycle$fault else NULL
        result <- tryCatch(helpers$runKitMember(helpers, owned, application_root, case$library, corpus, journey, fault, api),
            error = function(error) list(status = "driver_failed", error = conditionMessage(error)))
        # Counters are read while the slice may still exist; the monitor stops it on a guard.
        counters <- tryCatch(list(peak = api("workflowMemoryCounter")(file.path(group, "memory.peak")),
            swap = api("workflowMemoryCounter")(file.path(group, "memory.swap.peak")),
            accounting = api("workflowCgroupAccounting")(group)), error = function(error) NULL)
        file.create(stop_file)
        Sys.sleep(0.5)
        api("workflowSystemd")(c("stop", config$slice), required = FALSE)
        Sys.sleep(0.5)
        slice_state <- sliceState(config$slice)
        active_slice <- NULL
        trace <- utils::read.delim(trace_path, stringsAsFactors = FALSE)
        trace <- trace[is.finite(trace$uptime), , drop = FALSE]
        phases <- utils::read.delim(file.path(owned, "journey-phases.tsv"), stringsAsFactors = FALSE)
        guard_path <- file.path(application_root, "application-guard-stop.rds")
        guard <- if (file.exists(guard_path)) readRDS(guard_path) else NULL
        sampler_guard <- file.path(owned, "guard-stop.txt")
        outer_guard <- if (file.exists(sampler_guard)) readLines(sampler_guard, n = 1L) else NULL
        peak <- if (!is.null(counters)) counters$peak else max(trace$peak, na.rm = TRUE)
        swap <- if (!is.null(counters)) counters$swap else max(trace$swap_peak, na.rm = TRUE)
        head_unchanged <- NA
        if (case$expect_failure && !is.null(result$store) && !is.null(result$populations[[fault$after_step]])) {
            head <- helpers$kitReadHead(result$store, api)
            head_unchanged <- identical(head$head_sha256, result$populations[[fault$after_step]]$head_sha256)
        }
        operation <- c(list(status = if (!is.null(outer_guard)) outer_guard else result$status,
            config = config, binding = case$binding, member = case,
            elapsed_seconds = (result$teardown_clock %||% clock()) - operation_started,
            peak_bytes = peak, peak_source = if (!is.null(counters)) "cgroup" else "trace",
            swap_bytes = swap, accounting = counters$accounting,
            sampled_anon_bytes = max(trace$anon, na.rm = TRUE),
            max_sample_interval = if (nrow(trace) > 1L) max(diff(trace$uptime)) else NA_real_,
            samples = nrow(trace), guard = guard, outer_guard = outer_guard,
            head_unchanged = head_unchanged, populations_complete = isTRUE(result$populations_complete),
            cleanup = list(slice_stopped = identical(slice_state, "inactive") || !dir.exists(group)),
            boundary = "managed gui launch (controller, monitor, workers) through journey teardown; coordinator, browser and verifier outside the slice",
            clock = "shared kernel uptime from the launcher stamp",
            host_at_start = host, host_at_end = hostReady()),
            helpers$kitPhaseMeasures(phases, trace))
        saveRDS(operation, file.path(owned, "operation-receipt.rds"))
        receipt$cases[[case$id]] <- list(passed = operation$status %in% c("complete", "journey_failed"),
            operation = operation, result_status = result$status, error = result$error)
        saveRDS(receipt, receipt_path)
        verify_dir <- file.path(owned, "verify")
        dir.create(verify_dir)
        log <- file.path(owned, "verification-launch.log")
        verify_expression <- paste(
            'a <- commandArgs(TRUE); .libPaths(c(a[[1L]], .libPaths())); ns <- loadNamespace("MultiScholaR");',
            'api <- function(name) get(name, ns); env <- new.env(parent = globalenv());',
            'sys.source(file.path(a[[3L]], "tools/profiling/workflow_kit.R"), env);',
            'env$verifyKitMember(a[[2L]], a[[3L]], a[[4L]], api); cat("VERIFIED\\n")')
        status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", "-e", shQuote(verify_expression),
            shQuote(c(case$library, verify_dir, snapshot, case$selector))), stdout = log, stderr = log)
        verification_path <- file.path(owned, "verification.rds")
        verification <- if (file.exists(verification_path)) readRDS(verification_path) else NULL
        receipt$cases[[case$id]]$verifier_exit <- status
        receipt$cases[[case$id]]$verifier_completion <- list(accepted = identical(status, 0L) &&
            identical(verification$status, "complete"), reason = if (!identical(status, 0L)) paste(readLines(log, warn = FALSE), collapse = "\n"))
        saveRDS(receipt, receipt_path)
        cat("VERIFY", case$id, if (isTRUE(receipt$cases[[case$id]]$verifier_completion$accepted)) "passed" else "FAILED", "\n")
        cat("DONE", case$id, operation$status, "\n")
    }
    receipt$decision <- if (startsWith(stage, "kit-lifecycle-")) helpers$kitLifecycleDecision(receipt, cases, bindings) else
        list(status = "recorded", promotion_authorized = FALSE, note = "journey pair decisions are defined by KIT-005")
    saveRDS(receipt$decision, file.path(stage_root, "decision.rds"))
    receipt$status <- "complete"
}, error = function(error) {
    if (!is.null(sampler)) file.create(sampler)
    if (!is.null(active_slice)) api("workflowSystemd")(c("stop", active_slice), required = FALSE)
    receipt$status <<- if (identical(receipt$status, "setup_refused")) "setup_refused" else "failed"
    receipt$error <<- conditionMessage(error)
})
receipt$finished <- Sys.time()
saveRDS(receipt, receipt_path)
if (receipt$status != "complete") stop(receipt$error)
cat("DECISION", receipt$decision$status, receipt$decision$failing_rule %||% "", "\n")
