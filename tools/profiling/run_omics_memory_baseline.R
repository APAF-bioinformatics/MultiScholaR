#!/usr/bin/env Rscript

omicsMemoryScriptPath <- function() {
    file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(file_arg)) {
        return(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE))
    }
    source_paths <- unlist(lapply(sys.frames(), \(frame) {
        tryCatch(frame$ofile, error = \(...) NULL)
    }), use.names = FALSE)
    source_paths <- source_paths[grepl(
        "run_omics_memory_baseline[.]R$",
        source_paths
    )]
    if (length(source_paths)) {
        return(normalizePath(tail(source_paths, 1L), mustWork = FALSE))
    }
    candidates <- c(
        file.path(getwd(), "tools", "profiling", "run_omics_memory_baseline.R"),
        file.path(
            getwd(),
            "..",
            "..",
            "tools",
            "profiling",
            "run_omics_memory_baseline.R"
        )
    )
    existing <- candidates[file.exists(candidates)]
    if (!length(existing)) {
        stop("Cannot locate run_omics_memory_baseline.R", call. = FALSE)
    }
    normalizePath(existing[[1L]], mustWork = TRUE)
}

omicsMemoryRepoRoot <- function() {
    normalizePath(
        file.path(dirname(omicsMemoryScriptPath()), "..", ".."),
        mustWork = TRUE
    )
}

.OMICS_MEMORY_REPO_ROOT <- omicsMemoryRepoRoot()
source(
    file.path(
        .OMICS_MEMORY_REPO_ROOT,
        "tools",
        "profiling",
        "omics_resource_measurement.R"
    ),
    local = FALSE
)
source(
    file.path(
        .OMICS_MEMORY_REPO_ROOT,
        "tools",
        "profiling",
        "omics_workload_preparation.R"
    ),
    local = FALSE
)
source(
    file.path(
        .OMICS_MEMORY_REPO_ROOT,
        "tools",
        "profiling",
        "omics_workload_contract.R"
    ),
    local = FALSE
)

.OMICS_MEMORY_RESULT_SCHEMA <- "multischolar.omics_memory_baseline"
.OMICS_MEMORY_RESULT_VERSION <- "1.0.0"
.OMICS_MEMORY_WORKER_SCHEMA <- "multischolar.omics_memory_worker"
.OMICS_MEMORY_WORKER_VERSION <- "1.0.0"

omicsMemoryDefaultArgs <- function() {
    list(
        contract = NULL,
        adapter = NULL,
        output_dir = file.path("tests", "testthat", "_omics_memory_baseline"),
        diagnostics = TRUE,
        install_package = TRUE,
        prepare_spec = NULL,
        worker_spec = NULL,
        help = FALSE
    )
}

omicsMemoryParseBool <- function(value) {
    tolower(as.character(value)) %in% c("1", "true", "yes", "y")
}

omicsMemoryParseArgs <- function(argv) {
    args <- omicsMemoryDefaultArgs()
    index <- 1L
    while (index <= length(argv)) {
        token <- argv[[index]]
        if (identical(token, "--help")) {
            args$help <- TRUE
            index <- index + 1L
            next
        }
        if (!startsWith(token, "--")) {
            stop(sprintf("Unexpected positional argument: %s", token), call. = FALSE)
        }
        parts <- strsplit(sub("^--", "", token), "=", fixed = TRUE)[[1L]]
        key <- gsub("-", "_", parts[[1L]], fixed = TRUE)
        if (!key %in% names(args)) {
            stop(sprintf("Unknown option: %s", token), call. = FALSE)
        }
        if (length(parts) == 2L) {
            value <- parts[[2L]]
        } else {
            index <- index + 1L
            if (index > length(argv)) {
                stop(sprintf("Missing value for option: %s", token), call. = FALSE)
            }
            value <- argv[[index]]
        }
        if (key %in% c("diagnostics", "install_package")) {
            args[[key]] <- omicsMemoryParseBool(value)
        } else {
            args[[key]] <- value
        }
        index <- index + 1L
    }
    args
}

omicsMemoryUsage <- function() {
    cat(paste(
        "Usage: Rscript --vanilla tools/profiling/run_omics_memory_baseline.R [options]",
        "",
        "Options:",
        "  --contract <path>         Versioned workload contract",
        "  --adapter <path>          Omic-specific workload adapter",
        "  --output-dir <path>       Untracked evidence directory",
        "  --diagnostics <bool>      Run one separate allocation diagnostic",
        "  --install-package <bool>  Install the package once before measurement",
        "  --help",
        sep = "\n"
    ))
}

omicsMemoryResolvePath <- function(path, must_work = FALSE) {
    if (is.null(path)) {
        return(NULL)
    }
    resolved <- path.expand(path)
    if (!grepl("^(/|[A-Za-z]:[/\\\\])", resolved)) {
        resolved <- file.path(.OMICS_MEMORY_REPO_ROOT, resolved)
    }
    normalizePath(resolved, mustWork = must_work)
}

omicsMemoryEnvironment <- function(contract) {
    packages <- c("MultiScholaR", "processx", "ps", "jsonlite", "digest")
    versions <- lapply(packages, \(package) {
        if (!requireNamespace(package, quietly = TRUE)) {
            return(NULL)
        }
        as.character(utils::packageVersion(package))
    })
    names(versions) <- packages
    if (is.null(versions$MultiScholaR)) {
        versions$MultiScholaR <- read.dcf(
            file.path(.OMICS_MEMORY_REPO_ROOT, "DESCRIPTION"),
            fields = "Version"
        )[[1L]]
    }
    os_info <- Sys.info()
    public_os_fields <- intersect(
        c("sysname", "release", "version", "machine"),
        names(os_info)
    )
    list(
        captured_at = omicsResourceUtcNow(),
        code_revision = omicsResourceCodeRevision(.OMICS_MEMORY_REPO_ROOT),
        r_version = R.version.string,
        platform = R.version$platform,
        os = as.list(os_info[public_os_fields]),
        locale = contract$execution$locale,
        timezone = contract$execution$timezone,
        threads = as.integer(contract$execution$threads),
        rng_kind = as.list(unlist(contract$rng$kind, use.names = FALSE)),
        rng_seed = as.integer(contract$rng$seed),
        package_versions = versions,
        parent_thread_environment = omicsResourceThreadEnvironment()
    )
}

omicsMemoryDiskCategories <- function() {
    list(
        list(category = "diagnostics", pattern = "(^|/)diagnostics(/|$)"),
        list(
            category = "staging_snapshot",
            pattern = "(^|/)(staging|snapshots)(/|$)"
        ),
        list(category = "duckdb_spill", pattern = "(^|/)duckdb-spill(/|$)"),
        list(category = "committed", pattern = "(^|/)committed(/|$)"),
        list(category = "final", pattern = "(^|/)final(/|$)"),
        list(category = "harness", pattern = ".*")
    )
}

omicsMemoryValidateAdapterRun <- function(result) {
    omicsWorkloadRequireNames(
        result,
        c(
            "status",
            "workflow_evidence",
            "query_evidence",
            "session_evidence",
            "report_evidence",
            "retained"
        ),
        "adapter run result"
    )
    if (!identical(result$status, "passed")) {
        omicsWorkloadAbort("adapter run status must be 'passed'")
    }
    evidence_fields <- c(
        "workflow_evidence",
        "query_evidence",
        "session_evidence",
        "report_evidence"
    )
    if (any(!vapply(result[evidence_fields], is.list, logical(1)))) {
        omicsWorkloadAbort("all adapter evidence entries must be lists")
    }
    if (!isTRUE(result$workflow_evidence$truth_valid)) {
        omicsWorkloadAbort("adapter workflow evidence must prove truth_valid")
    }
    if (!length(result$query_evidence)) {
        omicsWorkloadAbort("adapter run must return bounded query evidence")
    }
    for (query in result$query_evidence) {
        required <- c("query_id", "rows", "p95_seconds")
        valid <- all(required %in% names(query)) &&
            is.character(query$query_id) && length(query$query_id) == 1L &&
            nzchar(query$query_id) &&
            is.finite(as.numeric(query$rows)) && as.numeric(query$rows) >= 0 &&
            is.finite(as.numeric(query$p95_seconds)) &&
            as.numeric(query$p95_seconds) >= 0
        if (!valid) {
            omicsWorkloadAbort("query evidence is missing finite required metrics")
        }
    }
    session_required <- c("status", "resource_count")
    report_required <- c("status", "file_count")
    session_valid <- all(session_required %in% names(result$session_evidence)) &&
        is.finite(as.numeric(result$session_evidence$resource_count)) &&
        as.numeric(result$session_evidence$resource_count) >= 0
    report_valid <- all(report_required %in% names(result$report_evidence)) &&
        is.finite(as.numeric(result$report_evidence$file_count)) &&
        as.numeric(result$report_evidence$file_count) >= 0
    if (!session_valid || !report_valid) {
        omicsWorkloadAbort("session and report evidence require finite counts")
    }
    invisible(result)
}

omicsMemoryDeterministicEvidence <- function(value) {
    if (!is.list(value)) {
        return(value)
    }
    if (!is.null(names(value))) {
        volatile <- grepl(
            "(^|_)(seconds|milliseconds|rss_bytes|vms_bytes|timings?)$",
            names(value)
        )
        value <- value[!volatile]
    }
    lapply(value, omicsMemoryDeterministicEvidence)
}

omicsMemoryWorker <- function(spec_path) {
    spec <- jsonlite::read_json(spec_path, simplifyVector = FALSE)
    contract <- omicsWorkloadReadContract(spec$contract_path)
    adapter <- omicsWorkloadLoadAdapter(spec$adapter_path, contract)
    rng_kind <- unlist(contract$rng$kind, use.names = FALSE)
    do.call(RNGkind, as.list(rng_kind))
    set.seed(as.integer(contract$rng$seed))
    Sys.setenv(
        TZ = contract$execution$timezone,
        OMP_NUM_THREADS = as.character(contract$execution$threads),
        OPENBLAS_NUM_THREADS = as.character(contract$execution$threads),
        MKL_NUM_THREADS = as.character(contract$execution$threads),
        ARROW_NUM_THREADS = as.character(contract$execution$threads),
        DUCKDB_THREADS = as.character(contract$execution$threads)
    )
    try(Sys.setlocale("LC_ALL", contract$execution$locale), silent = TRUE)
    if (!is.null(spec$package_library)) {
        loadNamespace("MultiScholaR", lib.loc = spec$package_library)
    }

    allocation_path <- file.path(spec$run_dir, "diagnostics", "rprofmem.out")
    if (isTRUE(spec$diagnostics_only)) {
        dir.create(dirname(allocation_path), recursive = TRUE, showWarnings = FALSE)
        utils::Rprofmem(allocation_path)
        on.exit(utils::Rprofmem(NULL), add = TRUE)
    }
    adapter_result <- adapter$run(list(
        contract = contract,
        payload_path = spec$payload_path,
        truth_path = spec$truth_path,
        run_dir = spec$run_dir,
        package_library = spec$package_library,
        diagnostics_only = isTRUE(spec$diagnostics_only)
    ))
    omicsMemoryValidateAdapterRun(adapter_result)
    observed_payload <- omicsWorkloadFileDigest(spec$payload_path)
    observed_truth <- omicsWorkloadFileDigest(spec$truth_path)
    binding_valid <- identical(
        observed_payload,
        contract$expected_digests$payload_sha256
    ) && identical(
        observed_truth,
        contract$expected_digests$truth_sha256
    ) && identical(
        spec$binding$binding_sha256,
        omicsWorkloadDigest(spec$binding[names(spec$binding) != "binding_sha256"])
    )
    if (!binding_valid) {
        omicsWorkloadAbort(
            "worker evidence binding does not match its inputs",
            class = "omics_workload_binding_error"
        )
    }
    if (isTRUE(spec$diagnostics_only)) {
        utils::Rprofmem(NULL)
        allocation <- omicsResourceAllocationMetrics(
            allocation_path,
            omicsWorkloadFileDigest
        )
    } else {
        allocation <- list(
            available = FALSE,
            reason = "allocation diagnostics run in a separate child process"
        )
    }
    result <- list(
        schema = .OMICS_MEMORY_WORKER_SCHEMA,
        schema_version = .OMICS_MEMORY_WORKER_VERSION,
        status = "passed",
        workload_id = contract$workload_id,
        capability_id = contract$capability$capability_id,
        binding_sha256 = spec$binding$binding_sha256,
        payload_sha256 = observed_payload,
        truth_sha256 = observed_truth,
        workflow_evidence = adapter_result$workflow_evidence,
        query_evidence = adapter_result$query_evidence,
        session_evidence = adapter_result$session_evidence,
        report_evidence = adapter_result$report_evidence,
        allocation_diagnostics = allocation,
        thread_environment = omicsResourceThreadEnvironment(),
        native_resources = omicsResourceObservedNativeMetrics(),
        retained_rss_bytes = omicsResourceProcSelfRss()
    )
    result$workflow_sha256 <- omicsWorkloadDigest(omicsMemoryDeterministicEvidence(list(
        workflow_evidence = result$workflow_evidence,
        query_evidence = result$query_evidence,
        session_evidence = result$session_evidence,
        report_evidence = result$report_evidence
    )))
    omicsResourceWriteJson(result, spec$result_path)
    writeLines("ready", file.path(spec$run_dir, "retention-ready"))
    Sys.sleep(as.numeric(contract$execution$retention_window_ms) / 1000)
    invisible(result)
}

omicsMemoryWorkerSafely <- function(spec_path) {
    spec <- jsonlite::read_json(spec_path, simplifyVector = FALSE)
    tryCatch(
        omicsMemoryWorker(spec_path),
        error = \(condition) {
            result <- list(
                schema = .OMICS_MEMORY_WORKER_SCHEMA,
                schema_version = .OMICS_MEMORY_WORKER_VERSION,
                status = "failed",
                reason = conditionMessage(condition),
                condition_class = class(condition)[[1L]]
            )
            omicsResourceWriteJson(result, spec$result_path)
            invisible(result)
        }
    )
}

omicsMemoryValidateWorker <- function(worker, binding) {
    valid <- is.list(worker) &&
        identical(worker$schema, .OMICS_MEMORY_WORKER_SCHEMA) &&
        identical(worker$schema_version, .OMICS_MEMORY_WORKER_VERSION) &&
        identical(worker$status, "passed") &&
        identical(worker$binding_sha256, binding$binding_sha256) &&
        identical(worker$payload_sha256, binding$payload_sha256) &&
        identical(worker$truth_sha256, binding$truth_sha256)
    if (!valid) {
        return(list(valid = FALSE, reason = worker$reason %||% "invalid worker result"))
    }
    required_metrics <- c("retained_rss_bytes")
    finite <- all(vapply(worker[required_metrics], \(value) {
        length(value) == 1L && is.finite(as.numeric(value))
    }, logical(1)))
    list(
        valid = finite,
        reason = if (finite) NULL else "worker returned non-finite metrics"
    )
}

omicsMemoryUndeclaredFiles <- function(run_dir) {
    paths <- list.files(
        run_dir,
        recursive = TRUE,
        full.names = FALSE,
        all.files = TRUE,
        no.. = TRUE
    )
    top_level <- sub(paste0("\\", .Platform$file.sep, ".*$"), "", paths)
    allowed <- c(
        "committed",
        "diagnostics",
        "duckdb-spill",
        "final",
        "harness",
        "retention-ready",
        "snapshots",
        "staging"
    )
    paths[!top_level %in% allowed]
}

omicsMemoryRunOne <- function(
    contract,
    contract_path,
    adapter_path,
    prepared,
    output_dir,
    environment,
    repetition,
    package_library,
    diagnostics_only = FALSE
) {
    suffix <- if (diagnostics_only) "diagnostic" else sprintf("r%02d", repetition)
    run_dir <- file.path(output_dir, "runs", suffix)
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
    binding <- omicsWorkloadEvidenceBinding(
        contract = contract,
        adapter_path = adapter_path,
        payload_path = prepared$payload_path,
        truth_path = prepared$truth_path,
        code_revision = environment$code_revision,
        environment = environment[c(
            "r_version",
            "platform",
            "locale",
            "timezone",
            "threads",
            "rng_kind",
            "rng_seed",
            "package_versions"
        )]
    )
    result_path <- file.path(run_dir, "final", "worker-result.json")
    spec_path <- tempfile("omics-memory-worker-", fileext = ".json")
    on.exit(unlink(spec_path, force = TRUE), add = TRUE)
    spec <- list(
        contract_path = contract_path,
        adapter_path = adapter_path,
        payload_path = prepared$payload_path,
        truth_path = prepared$truth_path,
        run_dir = normalizePath(run_dir, mustWork = TRUE),
        result_path = result_path,
        package_library = package_library,
        binding = binding,
        diagnostics_only = diagnostics_only
    )
    omicsResourceWriteJson(spec, spec_path)
    execution <- list(
        sampling_interval_ms = as.integer(contract$execution$sampling_interval_ms)
    )
    measured <- omicsResourceMeasureSubprocess(
        command = file.path(R.home("bin"), "Rscript"),
        command_args = c(
            "--vanilla",
            omicsMemoryScriptPath(),
            "--worker-spec",
            spec_path
        ),
        run_dir = run_dir,
        execution = execution,
        categories = omicsMemoryDiskCategories(),
        repo_root = .OMICS_MEMORY_REPO_ROOT,
        env = c(
            OMP_NUM_THREADS = as.character(contract$execution$threads),
            OPENBLAS_NUM_THREADS = as.character(contract$execution$threads),
            MKL_NUM_THREADS = as.character(contract$execution$threads),
            ARROW_NUM_THREADS = as.character(contract$execution$threads),
            DUCKDB_THREADS = as.character(contract$execution$threads),
            TZ = contract$execution$timezone
        ),
        timeout_ms = as.numeric(contract$execution$timeout_ms),
        maximum_retained_samples = as.integer(
            contract$execution$maximum_retained_samples
        ),
        require_retention_marker = TRUE
    )
    worker <- if (file.exists(result_path)) {
        jsonlite::read_json(result_path, simplifyVector = FALSE)
    } else {
        list(status = "missing", reason = "worker result was not written")
    }
    worker_validation <- omicsMemoryValidateWorker(worker, binding)
    undeclared <- omicsMemoryUndeclaredFiles(run_dir)
    metrics_finite <- all(vapply(
        measured$metrics[c(
            "elapsed_seconds",
            "peak_tree_rss_bytes",
            "retained_tree_rss_bytes",
            "peak_disk_bytes",
            "peak_artifact_disk_bytes",
            "final_disk_bytes",
            "final_file_count",
            "sample_count"
        )],
        \(value) length(value) == 1L && is.finite(as.numeric(value)),
        logical(1)
    ))
    passed <- identical(measured$status, "passed") &&
        isTRUE(worker_validation$valid) && metrics_finite && !length(undeclared)
    measured$status <- if (passed) "passed" else "failed"
    measured$repetition <- repetition
    measured$measurement_class <- if (diagnostics_only) {
        "allocation_copy_diagnostic"
    } else {
        "release_performance"
    }
    measured$cache_state <- contract$execution$cache_state
    measured$binding <- binding
    measured$preparation <- list(
        outside_measured_interval = TRUE,
        fresh_process = TRUE,
        metadata = prepared$metadata,
        payload_bytes = as.numeric(file.info(prepared$payload_path)$size),
        truth_bytes = as.numeric(file.info(prepared$truth_path)$size)
    )
    measured$metrics$committed_input_bytes <-
        measured$preparation$payload_bytes + measured$preparation$truth_bytes
    measured$metrics$committed_input_file_count <- 2L
    measured$worker <- worker
    measured$validation <- list(
        worker = worker_validation,
        finite_release_metrics = metrics_finite,
        undeclared_files = as.list(undeclared)
    )
    measured
}

omicsMemorySummarize <- function(runs) {
    passed <- vapply(runs, \(run) identical(run$status, "passed"), logical(1))
    list(
        requested = length(runs),
        completed = sum(passed),
        failed = sum(!passed),
        peak_tree_rss_bytes = as.list(vapply(
            runs,
            \(run) as.numeric(run$metrics$peak_tree_rss_bytes),
            numeric(1)
        )),
        retained_tree_rss_bytes = as.list(vapply(
            runs,
            \(run) as.numeric(run$metrics$retained_tree_rss_bytes),
            numeric(1)
        )),
        elapsed_seconds = as.list(vapply(
            runs,
            \(run) as.numeric(run$metrics$elapsed_seconds),
            numeric(1)
        ))
    )
}

omicsMemoryMain <- function(argv = commandArgs(trailingOnly = TRUE)) {
    args <- omicsMemoryParseArgs(argv)
    if (isTRUE(args$help)) {
        omicsMemoryUsage()
        return(invisible(0L))
    }
    if (!is.null(args$worker_spec)) {
        result <- omicsMemoryWorkerSafely(args$worker_spec)
        return(invisible(if (identical(result$status, "passed")) 0L else 1L))
    }
    if (!is.null(args$prepare_spec)) {
        result <- omicsMemoryPreparationWorkerSafely(args$prepare_spec)
        return(invisible(if (identical(result$status, "passed")) 0L else 1L))
    }
    if (is.null(args$contract) || is.null(args$adapter)) {
        stop("--contract and --adapter are required", call. = FALSE)
    }
    contract_path <- omicsMemoryResolvePath(args$contract, must_work = TRUE)
    adapter_path <- omicsMemoryResolvePath(args$adapter, must_work = TRUE)
    output_dir <- omicsMemoryResolvePath(args$output_dir)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    contract <- omicsWorkloadReadContract(contract_path)
    omicsWorkloadLoadAdapter(adapter_path, contract)
    environment <- omicsMemoryEnvironment(contract)
    preparations <- lapply(seq_len(2L), \(preparation_index) {
        omicsMemoryPrepareFresh(
            contract = contract,
            contract_path = contract_path,
            adapter_path = adapter_path,
            output_dir = output_dir,
            preparation_index = preparation_index
        )
    })
    preparation_deterministic <- identical(
        preparations[[1L]]$payload_sha256,
        preparations[[2L]]$payload_sha256
    ) && identical(
        preparations[[1L]]$truth_sha256,
        preparations[[2L]]$truth_sha256
    )
    if (!preparation_deterministic) {
        omicsWorkloadAbort(
            "fresh workload preparations produced different digests",
            class = "omics_workload_binding_error"
        )
    }
    package_library <- if (isTRUE(args$install_package)) {
        omicsResourceInstallPackage(.OMICS_MEMORY_REPO_ROOT, output_dir)
    } else {
        NULL
    }
    if (!is.null(package_library)) {
        on.exit(unlink(package_library, recursive = TRUE, force = TRUE), add = TRUE)
    }

    runs <- lapply(seq_len(as.integer(contract$execution$repetitions)), \(repetition) {
        omicsMemoryRunOne(
            contract = contract,
            contract_path = contract_path,
            adapter_path = adapter_path,
            prepared = preparations[[1L]],
            output_dir = output_dir,
            environment = environment,
            repetition = repetition,
            package_library = package_library
        )
    })
    diagnostics <- if (isTRUE(args$diagnostics)) {
        list(omicsMemoryRunOne(
            contract = contract,
            contract_path = contract_path,
            adapter_path = adapter_path,
            prepared = preparations[[1L]],
            output_dir = output_dir,
            environment = environment,
            repetition = 1L,
            package_library = package_library,
            diagnostics_only = TRUE
        ))
    } else {
        list()
    }
    payload_digests <- vapply(
        runs,
        \(run) run$binding$payload_sha256,
        character(1)
    )
    truth_digests <- vapply(
        runs,
        \(run) run$binding$truth_sha256,
        character(1)
    )
    workflow_digests <- vapply(runs, \(run) {
        run$worker$workflow_sha256 %||% NA_character_
    }, character(1))
    determinism <- list(
        fresh_preparation = preparation_deterministic,
        payload = length(unique(payload_digests)) == 1L,
        truth = length(unique(truth_digests)) == 1L,
        workflow = !anyNA(workflow_digests) &&
            length(unique(workflow_digests)) == 1L,
        binding = length(unique(vapply(
            runs,
            \(run) run$binding$binding_sha256,
            character(1)
        ))) == 1L
    )
    summary <- omicsMemorySummarize(runs)
    failed <- summary$failed > 0L || !all(unlist(determinism, use.names = FALSE))
    result <- list(
        schema = .OMICS_MEMORY_RESULT_SCHEMA,
        schema_version = .OMICS_MEMORY_RESULT_VERSION,
        status = if (failed) "failed" else "passed",
        generated_at = omicsResourceUtcNow(),
        workload_id = contract$workload_id,
        capability = contract$capability,
        contract_sha256 = omicsWorkloadDigest(contract),
        adapter = contract$adapter,
        environment = environment,
        measurement_contract = list(
            process_isolation = "one fresh measured child at a time",
            process_tree = "recursive RSS and VMS",
            retained_marker = "retention-ready",
            cache_state = contract$execution$cache_state,
            cache_claim = "uncontrolled unless the operating system cache was dropped",
            allocation_diagnostics = "separate child only",
            generated_inputs = "prepared before measured intervals",
            maximum_retained_samples = contract$execution$maximum_retained_samples
        ),
        preparation = list(
            fresh_processes = 2L,
            deterministic = preparation_deterministic,
            payload_sha256 = preparations[[1L]]$payload_sha256,
            truth_sha256 = preparations[[1L]]$truth_sha256,
            metadata = preparations[[1L]]$metadata
        ),
        determinism = determinism,
        runs = runs,
        diagnostics = diagnostics,
        summary = summary
    )
    result_path <- file.path(output_dir, "baseline-result.json")
    omicsResourceWriteJson(result, result_path)
    cat(sprintf(
        "Omics memory baseline status=%s runs=%d output=%s\n",
        result$status,
        length(runs),
        result_path
    ))
    invisible(if (failed) 1L else 0L)
}

`%||%` <- function(left, right) {
    if (is.null(left)) right else left
}

.OMICS_MEMORY_DIRECT <- any(
    normalizePath(
        sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)),
        mustWork = FALSE
    ) == omicsMemoryScriptPath()
)
if (isTRUE(.OMICS_MEMORY_DIRECT)) {
    status <- omicsMemoryMain()
    if (is.numeric(status) && status != 0L) {
        quit(status = status)
    }
}
