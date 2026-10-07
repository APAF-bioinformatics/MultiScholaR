# Explicit v2 qualification coordinator; scientific work runs only in installed workers.
args <- commandArgs(trailingOnly = TRUE)
# CORE stages accept a trailing --validate-only: bind and check everything, launch nothing.
validate_only <- length(args) == 4L && identical(args[[4L]], "--validate-only")
if (validate_only) args <- args[1:3]
core <- length(args) == 3L &&
    grepl("^core-(lifecycle-(complete|fail)|0(0[6-9]|1[0-4])|011-overflow)$", args[[1L]])
stopifnot(!validate_only || core)
stopifnot(length(args) == 3L, core || args[[1L]] %in%
    c("preflight-1", "preflight-2", "correction-3", "correction-4", "screen", "confirmation",
        "p14", "p14-diagnostic", "runtime-screen", "runtime-full",
        "runtime-repair-screen", "runtime-repair-screen-recheck", "runtime-repair-full",
        "runtime-staged-screen", "runtime-staged-full",
        "runtime-attribution-screen", "runtime-attribution-full", "runtime-batch-screen", "runtime-batch-full",
        "cache-screen", "cache-rollup", "cache-filters", "closure-filters",
        "structural-screen", "structural-rollup", "structural-primary", "structural-wide",
        "ordered-screen", "ordered-rollup", "ordered-primary", "ordered-wide"))
stage <- args[[1L]]
ordered <- startsWith(stage, "ordered-")
structural <- ordered || startsWith(stage, "structural-")
closure <- stage == "closure-filters"
cache_suite <- closure || startsWith(stage, "cache-")
attribution <- startsWith(stage, "runtime-attribution-")
batch <- startsWith(stage, "runtime-batch-")
runtime_staged <- stage %in% c("runtime-staged-screen", "runtime-staged-full")
runtime_repair <- stage %in% c("runtime-repair-screen", "runtime-repair-screen-recheck",
    "runtime-repair-full")
runtime_comparison <- structural || cache_suite || batch || attribution || runtime_staged || runtime_repair || stage %in% c("runtime-screen", "runtime-full")
profiling <- attribution || stage == "p14-diagnostic"
uptime_clock <- runtime_comparison || core
p14_stage <- stage %in% c("p14", "p14-diagnostic")
timed <- core || runtime_comparison || stage %in% c("screen", "confirmation", "p14", "p14-diagnostic")
member_seconds <- if (stage %in% c("confirmation", "p14", "runtime-full",
    "runtime-repair-full", "runtime-staged-full", "runtime-attribution-full", "runtime-batch-full",
    "cache-rollup", "cache-filters", "closure-filters", "structural-rollup", "structural-primary", "structural-wide",
    "ordered-rollup", "ordered-primary", "ordered-wide"))
    1800L else 300L
operationClock <- function() {
    if (p14_stage || uptime_clock) as.numeric(strsplit(readLines("/proc/uptime", n = 1L), " ")[[1L]][[1L]]) else
        proc.time()[["elapsed"]]
}
library_path <- normalizePath(args[[2L]])
root <- args[[3L]]
repo <- normalizePath(".")
.libPaths(c(library_path, .libPaths()))
loadNamespace("MultiScholaR")
api <- function(name) get(name, asNamespace("MultiScholaR"))
dir.create(root, recursive = TRUE, showWarnings = FALSE)
# A validate-only pass works in a private directory so it never claims the real stage.
stage_root <- if (validate_only) file.path(tempfile("core-validate-"), stage) else file.path(root, stage)
if (validate_only) dir.create(dirname(stage_root))
stopifnot(!file.exists(stage_root), dir.create(stage_root))
started <- Sys.time()
stage_seconds <- if (startsWith(stage, "correction-")) 300L else 600L
deadline <- started + stage_seconds
receipt_path <- file.path(stage_root, "stage.rds")
receipt <- list(stage = stage, status = "starting", started = started, cases = list(), diagnostic_only = profiling)
saveRDS(receipt, receipt_path)
snapshot <- file.path(stage_root, "source")
files <- c("tests/testthat/test-workflow-execution-adapters.R",
    "tools/profiling/workflow_perf106_repair.R",
    "tools/profiling/workflow_perf106.R",
    "tests/testthat/helper-workflow-checkpoint-fixtures.R",
    "tests/testthat/helper-workflow-bulk-rollup.R",
    "tests/testdata/workflow-execution/qualification.json",
    "tests/testdata/workflow-checkpoints/scenarios.json",
    "tickets/openai/verification/perf106-repair-protocol-v2.json")
if (runtime_comparison) files <- c(files, "tools/profiling/workflow_perf108.R",
    if (closure) "docs/qa/workflow-execution/perf-108/closeout/bindings.rds" else
    if (ordered) "docs/qa/workflow-execution/perf-108/ordered-groups/bindings.rds" else
        if (structural) "docs/qa/workflow-execution/perf-108/structural/bindings.rds" else
        if (cache_suite) "docs/qa/workflow-execution/perf-108/cache-qualification/bindings.rds" else
        if (batch) "docs/qa/workflow-execution/perf-108/runtime-batch-bindings.rds" else
        if (attribution) "docs/qa/workflow-execution/perf-108/runtime-attribution-bindings.rds" else
        if (runtime_staged) "docs/qa/workflow-execution/perf-108/runtime-staged-bindings.rds" else
        if (runtime_repair) "docs/qa/workflow-execution/perf-108/runtime-repair-bindings.rds" else
        "docs/qa/workflow-execution/perf-108/runtime-bindings.rds")
if (structural) files <- c(files, "tools/profiling/workflow_perf108_structural.R")
if (core) files <- c(files, "tools/profiling/workflow_perf108.R", "tools/profiling/workflow_core.R",
    "docs/qa/workflow-execution/core/bindings.json", "docs/qa/workflow-execution/core/protocol.md")
for (file in files) {
    target <- file.path(snapshot, file)
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    stopifnot(file.copy(file.path(repo, file), target))
}
receipt$source_hashes <- stats::setNames(vapply(file.path(snapshot, files),
    \(file) digest::digest(file = file, algo = "sha256"), character(1)), files)
helpers <- new.env(parent = baseenv())
sys.source(file.path(snapshot, "tools/profiling/workflow_perf106.R"), helpers)
if (runtime_comparison) sys.source(file.path(snapshot, "tools/profiling/workflow_perf108.R"), helpers)
if (structural) sys.source(file.path(snapshot, "tools/profiling/workflow_perf108_structural.R"), helpers)
if (core) sys.source(file.path(snapshot, "tools/profiling/workflow_core.R"), helpers)
if (core) {
    bindings_path <- file.path(snapshot, "docs/qa/workflow-execution/core/bindings.json")
    bindings <- helpers$coreBindings(bindings_path)
    cases <- helpers$coreCases(bindings, stage)
    member_seconds <- helpers$coreMemberSeconds(stage)
    frozen <- jsonlite::read_json(bindings_path)$harness
    receipt$input_bindings_sha256 <- digest::digest(file = bindings_path, algo = "sha256")
    receipt$build <- api("workflowExecutionBuildIdentity")()
    stopifnot(identical(normalizePath(root), normalizePath(bindings$root)),
        identical(digest::digest(file = file.path(snapshot, "docs/qa/workflow-execution/core/protocol.md"),
            algo = "sha256"), bindings$protocol_sha256),
        identical(unname(unlist(frozen$files)), vapply(file.path(repo, names(frozen$files)),
            \(file) digest::digest(file = file, algo = "sha256"), character(1), USE.NAMES = FALSE)),
        normalizePath(library_path) %in% normalizePath(vapply(cases, `[[`, character(1), "library")))
    for (case in cases) stopifnot(identical(digest::digest(file = file.path(case$input$store,
        "dia-state-head.rds"), algo = "sha256"), case$input$head_sha256), dir.exists(case$library))
    # Plan precondition (user 2026-10-08): 12 GiB startup headroom and no other managed unit.
    coreHostReady <- function() {
        units <- system2("systemctl", c("--user", "list-units", "--state=active", "--no-legend",
            "--plain", shQuote("multischolar*")), stdout = TRUE)
        list(available = api("workflowHostAvailable")(), units = units[nzchar(trimws(units))])
    }
    receipt$host_at_start <- coreHostReady()
    if (validate_only) {
        cat("VALIDATED", stage, length(cases), "case(s); host", format(receipt$host_at_start$available / 1024^3,
            digits = 3), "GiB available;", length(receipt$host_at_start$units), "active managed unit(s)\n")
        unlink(dirname(stage_root), recursive = TRUE)
        quit(save = "no", status = 0L)
    }
} else if (timed) {
    bindings_path <- if (runtime_comparison) file.path(snapshot,
        "docs/qa/workflow-execution/perf-108",
        if (closure) "closeout/bindings.rds" else if (ordered) "ordered-groups/bindings.rds" else
            if (structural) "structural/bindings.rds" else if (cache_suite) "cache-qualification/bindings.rds" else
            if (batch) "runtime-batch-bindings.rds" else if (attribution) "runtime-attribution-bindings.rds" else
            if (runtime_staged) "runtime-staged-bindings.rds" else
            if (runtime_repair) "runtime-repair-bindings.rds" else "runtime-bindings.rds") else
        file.path(root, "qualification-inputs.rds")
    bindings <- readRDS(bindings_path)
    stopifnot(identical(bindings$build, api("workflowExecutionBuildIdentity")()),
        identical(normalizePath(bindings$library), library_path))
    receipt$input_bindings_sha256 <- digest::digest(file = bindings_path, algo = "sha256")
    receipt$build <- bindings$build
    receipt$expected_bindings <- list()
    cases <- list()
    finalist <- NULL
    if (runtime_comparison) {
        stopifnot(identical(normalizePath(root), normalizePath(bindings$root)))
        if (closure) {
            stopifnot(isTRUE(bindings$retained_rollup_eligible),
                identical(bindings$worker_options, list()),
                identical(digest::digest(file = bindings$reused_rollup$path, algo = "sha256"),
                    bindings$reused_rollup$sha256))
        }
        if (structural && !stage %in% c("structural-screen", "ordered-screen")) stopifnot(isTRUE(readRDS(file.path(root,
            if (ordered) "ordered-screen/decision.rds" else "structural-screen/decision.rds"))$science_outputs_and_measurements_valid))
        if (stage == "cache-rollup") stopifnot(isTRUE(readRDS(file.path(root,
            "cache-screen/decision.rds"))$science_outputs_and_measurements_valid))
        if (stage == "cache-filters") stopifnot(identical(readRDS(file.path(root,
            "cache-rollup/decision.rds"))$status, "rollup_route_eligible"))
        if (attribution) stopifnot(identical(readRDS(file.path(root,
            "attribution-preflight.rds"))$status, "complete"))
        if (stage == "runtime-attribution-full") stopifnot(isTRUE(readRDS(file.path(root,
            "runtime-attribution-screen/decision.rds"))$science_outputs_and_measurements_valid))
        if (stage == "runtime-batch-full") stopifnot(isTRUE(readRDS(file.path(root,
            "runtime-batch-screen/decision.rds"))$science_outputs_and_measurements_valid))
        if (stage == "runtime-full") stopifnot(isTRUE(readRDS(file.path(root,
            "runtime-screen/decision-canonical.rds"))$science_outputs_and_measurements_valid))
        if (stage == "runtime-repair-full") stopifnot(isTRUE(readRDS(file.path(root,
            "runtime-repair-screen/decision-final.rds"))$science_outputs_and_measurements_valid))
        if (stage == "runtime-staged-full") stopifnot(isTRUE(readRDS(file.path(root,
            "runtime-staged-screen/decision.rds"))$science_outputs_and_measurements_valid))
        cases <- if (structural) helpers$perf108StructuralCases(bindings, stage) else
            if (cache_suite) helpers$perf108CacheCases(bindings, stage) else
            if (batch || attribution || runtime_repair || runtime_staged) helpers$perf108RepairCases(bindings, stage) else
            helpers$perf108Cases(bindings, stage)
        receipt$p14_builds <- lapply(bindings[if (cache_suite || structural)
            c("accepted", "current", "candidate") else c("control", "candidate")], `[[`, "build")
    } else if (p14_stage) {
        if (profiling) stopifnot(identical(readRDS(file.path(root, "profile-check.rds"))$status, "complete"))
        p14_path <- file.path(repo, "docs/qa/workflow-execution/perf-106/p14-regression-bindings.rds")
        p14 <- readRDS(p14_path)
        stopifnot(identical(digest::digest(file = p14$accepted_receipt, algo = "sha256"),
            p14$accepted_receipt_sha256), identical(digest::digest(file = file.path(p14$source_store,
                "dia-state-head.rds"), algo = "sha256"), p14$source_head_sha256))
        receipt$p14_input_bindings_sha256 <- digest::digest(file = p14_path, algo = "sha256")
        control <- readRDS(file.path(root, "correction-4/p14_control/application.rds"))
        receipt$p14_builds <- list(control = control$build, candidate = bindings$build)
        input <- list(store = p14$source_store, head_sha256 = p14$source_head_sha256,
            reference = p14$input_reference, parameters = list())
        roles <- if (profiling) c("control", "candidate") else c("control", "candidate", "candidate", "control")
        for (i in seq_along(roles)) {
            role <- roles[[i]]
            rep <- c(1L, 1L, 2L, 2L)[[i]]
            cases[[i]] <- list(id = paste(role, rep, sep = "-"), variant = "p14",
                role = role, control = role == "control", fault = "none", rep = rep,
                profile_worker = profiling,
                binding = digest::digest(list(input, receipt$p14_builds[[role]]), algo = "sha256"),
                input = input, verification_slot = 20L + i)
        }
    } else {
    if (stage == "confirmation") {
        decision_path <- file.path(root, "screen/decision.rds")
        decision <- readRDS(decision_path)
        stopifnot(identical(decision$status, "screen_finalist"),
            decision$finalist %in% c("physical", "selection"))
        finalist <- decision$finalist
        receipt$screen_decision_sha256 <- digest::digest(file = decision_path, algo = "sha256")
        receipt$finalist <- finalist
    }
    for (population in c("primary", "wide")) {
        input <- bindings[[if (stage == "screen") "screen" else "full"]][[population]]
        repetitions <- c(native = 0L, physical = 0L, selection = 0L)
        order <- if (stage == "screen")
            c("native", "physical", "selection", "selection", "physical", "native") else
            c("native", finalist, finalist, "native")
        for (variant in order) {
            repetitions[[variant]] <- repetitions[[variant]] + 1L
            binding <- digest::digest(list(input, bindings$build, variant,
                protocol = receipt$source_hashes[["tickets/openai/verification/perf106-repair-protocol-v2.json"]]),
                algo = "sha256")
            receipt$expected_bindings[[population]][[variant]] <- binding
            cases[[length(cases) + 1L]] <- list(
                id = paste(population, variant, repetitions[[variant]], sep = "-"),
                variant = variant, fault = "none", population = population,
                rep = repetitions[[variant]], binding = binding, input = input,
                verification_slot = length(cases) + 1L + if (stage == "confirmation") 12L else 0L)
        }
    }
    }
} else cases <- if (stage == "correction-4") list(
    list(id = "p14_control", variant = "p14", fault = "none", control = TRUE)
) else if (stage == "correction-3") list(
    list(id = "summary", variant = "physical", fault = "none", summary_heavy = TRUE)
) else if (stage == "preflight-1") list(
    list(id = "native", variant = "native", fault = "none"),
    list(id = "physical", variant = "physical", fault = "none"),
    list(id = "selection", variant = "selection", fault = "none"),
    list(id = "p14", variant = "p14", fault = "none"),
    list(id = "summary", variant = "physical", fault = "none", summary_heavy = TRUE)
) else list(
    list(id = "publication", variant = "physical", fault = "publication"),
    list(id = "worker_failure", variant = "selection", fault = "worker_failure"),
    list(id = "cancel", variant = "physical", fault = "cancel"),
    list(id = "worker_loss", variant = "native", fault = "worker_loss"),
    list(id = "monitor_loss", variant = "physical", fault = "monitor_loss"),
    list(id = "p14_control", variant = "p14", fault = "none", control = TRUE)
)
active_slice <- NULL
tryCatch({
    for (case in cases) {
        properties <- list(available = FALSE, reason = "service not yet observed")
        outer_guard <- NULL
        trace <- NULL
        if (timed) deadline <- Sys.time() + member_seconds
        if (Sys.time() >= deadline) stop("installed preflight ceiling exhausted")
        if (core) {
            host <- coreHostReady()
            if (host$available < 12 * 1024^3 || length(host$units)) {
                outer_guard <- "setup_refused"
                stop("CORE setup refused: ", format(host$available / 1024^3, digits = 3),
                    " GiB available, ", length(host$units), " active managed unit(s)")
            }
        }
        owned <- file.path(stage_root, case$id)
        stopifnot(dir.create(owned))
        member <- c(list(schema = 1L, programme = "perf106-repair-v2", preflight = !timed,
            operation = if (case$variant == "p14") "precursor_rollup" else "minimum_peptides_per_sample",
            deadline = if (timed) member_seconds else stage_seconds,
            parameters = if (timed) case$input$parameters else if (case$variant == "p14") list() else
                list(peptides_per_sample_cutoff = 2L, inclusion_list = character())), case)
        if (timed) {
            member$source_store <- case$input$store
            member$source_head_sha256 <- case$input$head_sha256
            member$input_entity_id <- case$input$reference$entity_id
        }
        saveRDS(member, file.path(owned, "timed-member.rds"))
        operation_started <- operationClock()
        if (core || runtime_comparison || isTRUE(case$control) || p14_stage) {
            # Load the frozen predecessor in a fresh launcher process; never
            # replace functions in either installed scientific namespace.
            control_library <- if (core || runtime_comparison) case$library else if (isTRUE(case$control)) file.path(dirname(dirname(library_path)),
                "browser-smoke.LGgIYq", "library") else library_path
            expression <- paste(
                'a <- commandArgs(TRUE); .libPaths(c(a[[1L]], .libPaths()));',
                'ns <- loadNamespace("MultiScholaR");',
                'stamp <- as.numeric(strsplit(readLines("/proc/uptime", n = 1L), " ")[[1L]][[1L]]);',
                'saveRDS(stamp, file.path(a[[2L]], "launcher-start.rds"));',
                'get("launchWorkflowApplication", ns)(a[[2L]], diagnostic = FALSE,',
                'mode = "adapters", test_root = a[[3L]])')
            status <- system2(file.path(R.home("bin"), "Rscript"),
                c("--vanilla", "-e", shQuote(expression),
                    shQuote(c(control_library, owned, snapshot))))
            stopifnot(identical(status, 0L))
            config <- readRDS(file.path(owned, "application.rds"))
            if (p14_stage || uptime_clock) {
                operation_started <- readRDS(file.path(owned, "launcher-start.rds"))
                deadline <- Sys.time() + member_seconds - (operationClock() - operation_started)
            }
        } else config <- api("launchWorkflowApplication")(owned, diagnostic = FALSE,
            mode = "adapters", test_root = snapshot)
        cat("START", case$id, config$controller_unit, "\n")
        active_slice <- config$slice
        if (core || runtime_comparison) stopifnot(identical(config$build, case$expected_build))
        if (p14_stage) stopifnot(identical(config$build, receipt$p14_builds[[case$role]]))
        group <- file.path("/sys/fs/cgroup", sub("^/", "", config$control_group))
        outer_guard <- NULL
        max_anon <- max_file <- 0
        trace <- NULL
        if (profiling) {
            trace <- file(file.path(owned, "memory-trace.tsv"), "wt")
            writeLines("uptime\tcurrent\tpeak\tanon\tfile\tkernel\tpagetables\tslab", trace)
        }
        repeat {
            if (timed) {
                stats <- strsplit(readLines(file.path(group, "memory.stat")), " ", fixed = TRUE)
                stats <- stats::setNames(vapply(stats, \(x) as.numeric(x[[2L]]), numeric(1)),
                    vapply(stats, `[[`, character(1), 1L))
                max_anon <- max(max_anon, stats[["anon"]])
                max_file <- max(max_file, stats[["file"]])
                if (profiling) {
                    writeLines(paste(c(operationClock(),
                        api("workflowMemoryCounter")(file.path(group, "memory.current")),
                        api("workflowMemoryCounter")(file.path(group, "memory.peak")),
                        stats[c("anon", "file", "kernel", "pagetables", "slab")]), collapse = "\t"), trace)
                    flush(trace)
                }
                if (api("workflowMemoryCounter")(file.path(group, "memory.peak")) >= 15.5 * 1024^3)
                    outer_guard <- "memory_guard"
                if (api("workflowHostAvailable")() < 2 * 1024^3) outer_guard <- "host_guard"
                if (Sys.time() >= deadline) outer_guard <- "time_guard"
                if (!is.null(outer_guard)) break
            }
            properties <- api("workflowServiceProperties")(config$controller_unit)
            if (!isTRUE(properties$available) || !properties$SubState %in% c("running", "start", "start-pre")) break
            if (Sys.time() >= deadline) {
                api("workflowSystemd")(c("stop", config$slice), required = FALSE)
                stop("installed preflight ceiling exhausted")
            }
            Sys.sleep(if (profiling) 0.01 else if (timed) 0.05 else 0.5)
        }
        if (!is.null(trace)) { close(trace); trace <- NULL }
        if (timed) {
            accounting <- api("workflowCgroupAccounting")(group)
            peak <- api("workflowMemoryCounter")(file.path(group, "memory.peak"))
            swap <- api("workflowMemoryCounter")(file.path(group, "memory.swap.peak"))
        }
        member_result <- file.path(owned, "member-result.rds")
        tests_path <- file.path(owned, "lifecycle-tests.rds")
        result <- if (file.exists(member_result)) readRDS(member_result) else NULL
        tests <- if (file.exists(tests_path)) readRDS(tests_path) else NULL
        passed <- !is.null(result) && !is.null(tests) && !any(tests$failed > 0 | tests$error)
        receipt$cases[[case$id]] <- list(passed = passed, result = result, properties = properties,
            assertions = if (!is.null(tests)) sum(tests$passed) else 0L)
        saveRDS(receipt, receipt_path)
        api("workflowSystemd")(c("stop", config$slice), required = FALSE)
        active_slice <- NULL
        if (timed) {
            operation <- list(status = if (!is.null(outer_guard)) outer_guard else
                    if (is.null(result)) "controller_failed" else result$status,
                config = config, binding = case$binding, member = member,
                diagnostic_only = profiling,
                elapsed_seconds = operationClock() - operation_started,
                peak_bytes = peak, swap_bytes = swap, accounting = accounting,
                sampled_anon_bytes = max_anon, sampled_file_bytes = max_file,
                lifecycle_passed = passed,
                boundary = "managed launch through publication and controller/worker/monitor teardown; separate verifier excluded",
                coordinator_excluded = "metadata-only launcher outside managed slice",
                clock = if (p14_stage || uptime_clock) "shared kernel uptime; namespace bootstrap excluded for both builds" else
                    "parent elapsed clock",
                execution = if (!is.null(result)) result$execution else NULL)
            operation$measurement_valid <- isTRUE(result$execution$measurement_valid) &&
                peak >= result$execution$peak_bytes &&
                operation$elapsed_seconds >= result$execution$elapsed_seconds &&
                isTRUE(result$execution$cleanup$worker_stopped) &&
                isTRUE(result$execution$cleanup$monitor_stopped)
            operation_path <- file.path(owned, "operation-receipt.rds")
            saveRDS(operation, operation_path)
            receipt$cases[[case$id]]$operation <- operation
            saveRDS(receipt, receipt_path)
            if ((operation$status == "complete" || core && case$expect_failure) && passed) {
                if (profiling) {
                    request_dir <- file.path(owned, "project/.workflow-entities/.workflow-export",
                        operation$execution$request_id)
                    profile_start <- readRDS(file.path(request_dir, "profile-start.rds"))
                    profile_finish <- readRDS(file.path(request_dir, "profile-finished.rds"))
                    stopifnot(profile_finish$finished > profile_start$start,
                        !isTRUE(profile_finish$error),
                        file.size(file.path(request_dir, "worker-stack.out")) > 0,
                        file.exists(file.path(request_dir, "worker-allocations.out")))
                }
                selector <- if (core) case$selector else if (closure) paste0("perf108:E", case$verification_slot) else if (ordered) paste0("perf108:W", switch(stage,
                    `ordered-screen` = "S", `ordered-rollup` = "R", `ordered-primary` = "P", `ordered-wide` = "W"),
                    case$verification_slot) else if (structural) paste0("perf108:Y", switch(stage,
                    `structural-screen` = "S", `structural-rollup` = "R", `structural-primary` = "P",
                    `structural-wide` = "W"), case$verification_slot) else
                    if (cache_suite) paste0("perf108:X", switch(stage,
                    `cache-screen` = "S", `cache-rollup` = "R", `cache-filters` = "F"), case$verification_slot) else
                    if (batch) paste0("perf108:H",
                    if (stage == "runtime-batch-full") "F1" else "S1") else if (attribution) paste0("perf108:A",
                    if (stage == "runtime-attribution-full") "F1" else "S1") else
                    if (runtime_staged) paste0("perf108:P",
                    if (stage == "runtime-staged-full") "F" else "S", case$rep) else
                    if (runtime_repair) paste0("perf108:R",
                    if (stage == "runtime-repair-full") "F" else if
                        (stage == "runtime-repair-screen-recheck") "R" else "S", case$rep) else
                    if (runtime_comparison) paste0("perf108:V", case$verification_slot) else if (profiling) paste0("perf106:V2-D", case$verification_slot - 20L) else
                    paste0("perf106:V2-Q", case$verification_slot)
                log <- file.path(owned, "verification-launch.log")
                status <- system2("bash", shQuote(c(
                    file.path(repo, "tools/profiling/run_workflow_checkpoint_proof.sh"),
                    "--service", "reader-check", library_path, dirname(root), repo, selector)),
                    stdout = log, stderr = log)
                receipt$cases[[case$id]]$verifier_exit <- status
                receipt$cases[[case$id]]$verifier_completion <- helpers$readPerf106Verification(
                    owned, dirname(root), status)
                saveRDS(receipt, receipt_path)
                cat("VERIFY", case$id, if (isTRUE(receipt$cases[[case$id]]$verifier_completion$accepted))
                    "passed" else paste("FAILED:", receipt$cases[[case$id]]$verifier_completion$reason), "\n")
                if (profiling && !isTRUE(receipt$cases[[case$id]]$verifier_completion$accepted))
                    stop(receipt$cases[[case$id]]$verifier_completion$reason)
            }
        }
        cat("DONE", case$id, if (passed) "operation checks passed" else "operation checks FAILED", "\n")
        if (profiling && !passed) stop(paste("diagnostic case failed:", case$id))
        if (!passed && !timed) stop(paste("installed case failed:", case$id))
    }
    if (core) {
        receipt$decision <- if (startsWith(stage, "core-lifecycle-")) helpers$coreLifecycleDecision(receipt, cases) else
            helpers$coreDecision(receipt, cases, sub("core-", "", stage, fixed = TRUE))
        saveRDS(receipt$decision, file.path(stage_root, "decision.rds"))
    } else if (structural) {
        receipt$decision <- helpers$perf108StructuralDecision(receipt, cases)
        saveRDS(receipt$decision, file.path(stage_root, "decision.rds"))
    } else if (cache_suite) {
        receipt$decision <- helpers$perf108CacheDecision(receipt, cases)
        saveRDS(receipt$decision, file.path(stage_root, "decision.rds"))
    } else if (attribution || batch) {
        checked <- receipt$cases[["p14-candidate-1"]]
        receipt$decision <- list(status = if (batch) "exploratory_single_member" else "diagnostic_only", promotion_authorized = FALSE,
            qualification_unchanged = TRUE,
            science_outputs_and_measurements_valid = isTRUE(checked$verifier_completion$accepted) &&
                isTRUE(checked$operation$measurement_valid))
        saveRDS(receipt$decision, file.path(stage_root, "decision.rds"))
    } else if (runtime_comparison) {
        receipt$decision <- if (runtime_repair || runtime_staged)
            helpers$perf108RepairDecision(receipt, cases, bindings) else
            helpers$perf108Decision(receipt, cases, bindings)
        saveRDS(receipt$decision, file.path(stage_root, "decision.rds"))
    } else if (profiling) {
        receipt$decision <- list(status = "diagnostic_only", promotion_authorized = FALSE,
            qualification_unchanged = TRUE)
        saveRDS(receipt$decision, file.path(stage_root, "decision.rds"))
    } else if (stage == "p14") {
        members <- list()
        valid <- TRUE
        ratios <- list()
        for (case in cases) {
            op <- receipt$cases[[case$id]]$operation
            checked <- helpers$readPerf106Verification(file.path(stage_root, case$id),
                dirname(root), receipt$cases[[case$id]]$verifier_exit)
            members[[case$id]] <- list(operation = op,
                science = isTRUE(checked$accepted))
        }
        for (rep in 1:2) {
            control <- members[[paste("control", rep, sep = "-")]]
            candidate <- members[[paste("candidate", rep, sep = "-")]]
            acceptable <- function(x) isTRUE(x$science) &&
                identical(x$operation$status, "complete") && isTRUE(x$operation$measurement_valid) &&
                x$operation$peak_bytes < 15.5 * 1024^3 && x$operation$elapsed_seconds <= 1800 &&
                identical(x$operation$swap_bytes, 0)
            paired <- acceptable(control) && acceptable(candidate)
            ratio <- if (paired) list(peak = candidate$operation$peak_bytes / control$operation$peak_bytes,
                elapsed = candidate$operation$elapsed_seconds / control$operation$elapsed_seconds) else NULL
            ratios[[as.character(rep)]] <- list(valid = paired, ratios = ratio)
            valid <- valid && paired && ratio$peak <= 1.05 && ratio$elapsed <= 1.2
        }
        receipt$decision <- list(status = if (valid) "p14_regression_passed" else "p14_regression_not_qualified",
            pairs = ratios, members = members, promotion_authorized = FALSE)
        saveRDS(receipt$decision, file.path(stage_root, "decision.rds"))
    } else if (timed) {
        populations <- list(primary = list(), wide = list())
        fields <- c("rows", "samples", "columns", "page_rows", "page_hashes",
            "design_sha256", "summary_sha256", "support_sha256", "required_outputs_complete")
        for (population in names(populations)) {
            oracle_id <- paste(population, "native", 1L, sep = "-")
            checked_oracle <- helpers$readPerf106Verification(file.path(stage_root, oracle_id),
                dirname(root), receipt$cases[[oracle_id]]$verifier_exit)
            oracle <- if (isTRUE(checked_oracle$accepted)) checked_oracle$value else NULL
            for (arm in if (stage == "screen") c("native", "physical", "selection") else
                c("native", finalist)) {
                populations[[population]][[arm]] <- lapply(1:2, function(rep) {
                    id <- paste(population, arm, rep, sep = "-")
                    op <- receipt$cases[[id]]$operation
                    checked <- helpers$readPerf106Verification(file.path(stage_root, id),
                        dirname(root), receipt$cases[[id]]$verifier_exit)
                    verification <- if (isTRUE(checked$accepted)) checked$value else NULL
                    exact <- !is.null(oracle) && identical(oracle$status, "complete") &&
                        !is.null(verification) && identical(verification$status, "complete") &&
                        identical(verification$verification[fields], oracle$verification[fields])
                    list(status = op$status, bindings = op$binding, rep = rep,
                        request_id = op$execution$request_id %||% id,
                        science = exact, outputs = isTRUE(verification$verification$required_outputs_complete),
                        measurement_valid = isTRUE(op$measurement_valid),
                        integrity = isTRUE(checked$accepted),
                        guard = isTRUE(op$peak_bytes < 15.5 * 1024^3) && identical(op$swap_bytes, 0),
                        admission = identical(op$status, "complete"),
                        peak_bytes = op$peak_bytes, elapsed_seconds = op$elapsed_seconds,
                        swap_bytes = op$swap_bytes,
                        spill_bytes = if (identical(verification$status, "complete")) 0 else NA_real_)
                })
            }
        }
        if (stage == "screen") {
            receipt$decision <- api("rankPerfRepairScreen")(populations, receipt$expected_bindings)
        } else {
            ratios <- list()
            valid <- TRUE
            for (population in names(populations)) for (rep in 1:2) {
                native <- populations[[population]]$native[[rep]]
                candidate <- populations[[population]][[finalist]][[rep]]
                acceptable <- function(x) identical(x$status, "complete") &&
                    all(vapply(c("science", "outputs", "measurement_valid", "integrity", "guard", "admission"),
                        \(key) isTRUE(x[[key]]), logical(1))) && x$elapsed_seconds <= 1800
                paired <- acceptable(native) && acceptable(candidate)
                ratio <- if (paired) list(peak = candidate$peak_bytes / native$peak_bytes,
                    elapsed = candidate$elapsed_seconds / native$elapsed_seconds) else NULL
                ratios[[paste(population, rep, sep = "-")]] <- list(valid = paired, ratios = ratio)
                valid <- valid && paired && ratio$peak <= 0.8 && ratio$elapsed <= 1.2
            }
            receipt$decision <- list(status = if (valid) "confirmation_passed" else "confirmation_not_qualified",
                finalist = finalist, pairs = ratios, promotion_authorized = FALSE)
        }
        receipt$decision$members <- populations
        saveRDS(receipt$decision, file.path(stage_root, "decision.rds"))
    }
    receipt$status <- "complete"
}, error = function(error) {
    if (exists("trace") && !is.null(trace)) close(trace)
    if (!is.null(active_slice)) api("workflowSystemd")(c("stop", active_slice), required = FALSE)
    receipt$status <<- "failed"
    receipt$error <<- conditionMessage(error)
    if (timed && exists("owned") && dir.exists(owned) &&
        !file.exists(file.path(owned, "operation-receipt.rds")))
        saveRDS(list(status = if (exists("outer_guard") && !is.null(outer_guard)) outer_guard else
                "instrumentation_failed", error = conditionMessage(error),
            member = member), file.path(owned, "operation-receipt.rds"))
})
receipt$finished <- Sys.time()
saveRDS(receipt, receipt_path)
if (receipt$status != "complete") stop(receipt$error)
