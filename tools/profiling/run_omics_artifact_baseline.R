#!/usr/bin/env Rscript

baselineScriptPath <- function() {
    file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(file_arg)) {
        return(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE))
    }
    source_path <- tryCatch(sys.frame(1L)$ofile, error = \(...) NULL)
    if (!is.null(source_path)) {
        return(normalizePath(source_path, mustWork = FALSE))
    }
    normalizePath(
        file.path(getwd(), "tools", "profiling", "run_omics_artifact_baseline.R"),
        mustWork = FALSE
    )
}

baselineRepoRoot <- function() {
    normalizePath(file.path(dirname(baselineScriptPath()), "..", ".."), mustWork = TRUE)
}

.BASELINE_REPO_ROOT <- baselineRepoRoot()
if (!exists("omicsParityReadManifest", mode = "function")) {
    source(
        file.path(.BASELINE_REPO_ROOT, "tests", "testthat", "helper-omics-parity.R"),
        local = FALSE
    )
}

baselineUtcNow <- function() {
    format(as.POSIXct(Sys.time(), tz = "UTC"), "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC")
}

baselineParseBool <- function(value) {
    tolower(as.character(value)) %in% c("1", "true", "yes", "y")
}

baselineThreadEnvironment <- function() {
    as.list(Sys.getenv(
        c(
            "OMP_NUM_THREADS",
            "OPENBLAS_NUM_THREADS",
            "MKL_NUM_THREADS",
            "ARROW_NUM_THREADS",
            "DUCKDB_THREADS"
        ),
        unset = NA_character_
    ))
}

baselineResolvePath <- function(path, must_work = FALSE) {
    resolved <- path.expand(path)
    if (!grepl("^(/|[A-Za-z]:[/\\\\])", resolved)) {
        resolved <- file.path(.BASELINE_REPO_ROOT, resolved)
    }
    normalizePath(resolved, mustWork = must_work)
}

baselineDefaultArgs <- function() {
    list(
        manifest = file.path("tests", "testdata", "omics-parity", "scenarios.json"),
        output_dir = file.path("tests", "testthat", "_omics_artifact_baseline"),
        mode = "fixture",
        scenario = character(),
        repetitions = NULL,
        diagnostics = TRUE,
        include_private = FALSE,
        worker_spec = NULL,
        help = FALSE
    )
}

baselineUsage <- function() {
    cat(paste(
        "Usage: Rscript --vanilla tools/profiling/run_omics_artifact_baseline.R [options]",
        "",
        "Options:",
        "  --manifest <path>          Scenario manifest",
        "  --output-dir <path>        Untracked result directory",
        "  --mode <fixture|scientific|all>",
        "  --scenario <id>            Repeatable fixture scenario selector",
        "  --repetitions <n>          Override manifest repetitions",
        "  --diagnostics <bool>        Run allocation/copy diagnostics separately",
        "  --include-private <bool>   Enable the opt-in private scenario",
        "  --help",
        sep = "\n"
    ))
}

baselineParseArgs <- function(argv) {
    args <- baselineDefaultArgs()
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
        if (identical(key, "scenario")) {
            args$scenario <- c(args$scenario, value)
        } else if (identical(key, "repetitions")) {
            args$repetitions <- as.integer(value)
        } else if (key %in% c("include_private", "diagnostics")) {
            args[[key]] <- baselineParseBool(value)
        } else {
            args[[key]] <- value
        }
        index <- index + 1L
    }
    args
}

baselineWriteJson <- function(value, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    jsonlite::write_json(
        value,
        path,
        auto_unbox = TRUE,
        pretty = TRUE,
        null = "null",
        na = "null",
        digits = 17
    )
    invisible(path)
}

baselineDisplayCommandArgs <- function(command_args, run_dir) {
    replacements <- c(
        "<repo-root>" = .BASELINE_REPO_ROOT,
        "<run-dir>" = normalizePath(run_dir, mustWork = FALSE)
    )
    vapply(command_args, \(argument) {
        displayed <- argument
        for (label in names(replacements)) {
            displayed <- gsub(replacements[[label]], label, displayed, fixed = TRUE)
        }
        displayed
    }, character(1), USE.NAMES = FALSE)
}

baselineDirMetrics <- function(path, categories) {
    files <- if (dir.exists(path)) {
        list.files(path, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
    } else {
        character()
    }
    files <- files[file.exists(files) & !dir.exists(files)]
    relative <- if (length(files)) {
        substring(files, nchar(normalizePath(path)) + 2L)
    } else {
        character()
    }
    sizes <- if (length(files)) as.numeric(file.info(files)$size) else numeric()
    category_names <- vapply(categories, `[[`, character(1), "category")
    bytes <- stats::setNames(as.list(rep(0, length(categories))), category_names)
    counts <- stats::setNames(as.list(rep(0L, length(categories))), category_names)
    assigned <- rep(FALSE, length(files))

    for (category in categories) {
        matches <- !assigned & grepl(category$pattern, relative, perl = TRUE)
        bytes[[category$category]] <- sum(sizes[matches], na.rm = TRUE)
        counts[[category$category]] <- sum(matches)
        assigned[matches] <- TRUE
    }
    list(
        total_bytes = sum(sizes, na.rm = TRUE),
        file_count = length(files),
        bytes = bytes,
        counts = counts
    )
}

baselineProcessTreeMetrics <- function(pid) {
    if (!requireNamespace("ps", quietly = TRUE)) {
        stop("Package 'ps' is required for process-tree RSS measurement", call. = FALSE)
    }
    root <- tryCatch(ps::ps_handle(pid), error = \(...) NULL)
    if (is.null(root) || !isTRUE(tryCatch(ps::ps_is_running(root), error = \(...) FALSE))) {
        return(list(rss_bytes = 0, vms_bytes = 0, process_count = 0L, pids = list()))
    }
    handles <- c(
        list(root),
        tryCatch(ps::ps_children(root, recursive = TRUE), error = \(...) list())
    )
    memory <- lapply(handles, \(handle) {
        tryCatch(ps::ps_memory_info(handle), error = \(...) c(rss = 0, vms = 0))
    })
    list(
        rss_bytes = sum(vapply(memory, \(info) as.numeric(info[["rss"]]), numeric(1))),
        vms_bytes = sum(vapply(memory, \(info) as.numeric(info[["vms"]]), numeric(1))),
        process_count = length(handles),
        pids = as.list(vapply(handles, \(handle) ps::ps_pid(handle), integer(1)))
    )
}

baselineMeasureSubprocess <- function(
    command,
    command_args,
    run_dir,
    execution,
    categories,
    env = character(),
    capture_output = TRUE
) {
    if (!requireNamespace("processx", quietly = TRUE)) {
        stop("Package 'processx' is required for fresh-process measurement", call. = FALSE)
    }
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
    stdout_path <- file.path(run_dir, "final", "stdout.log")
    stderr_path <- file.path(run_dir, "final", "stderr.log")
    dir.create(dirname(stdout_path), recursive = TRUE, showWarnings = FALSE)
    stdout_destination <- if (capture_output) stdout_path else nullfile()
    stderr_destination <- if (capture_output) stderr_path else nullfile()
    started_at <- baselineUtcNow()
    started <- proc.time()[["elapsed"]]
    process <- processx::process$new(
        command,
        command_args,
        stdout = stdout_destination,
        stderr = stderr_destination,
        env = env,
        wd = .BASELINE_REPO_ROOT,
        cleanup = TRUE,
        supervise = TRUE
    )
    samples <- list()
    peak_category_bytes <- stats::setNames(
        rep(0, length(categories)),
        vapply(categories, `[[`, character(1), "category")
    )
    retained <- NULL
    index <- 0L

    repeat {
        index <- index + 1L
        tree <- baselineProcessTreeMetrics(process$get_pid())
        disk <- baselineDirMetrics(run_dir, categories)
        for (category in names(peak_category_bytes)) {
            peak_category_bytes[[category]] <- max(
                peak_category_bytes[[category]],
                disk$bytes[[category]]
            )
        }
        sample <- list(
            elapsed_seconds = proc.time()[["elapsed"]] - started,
            tree_rss_bytes = tree$rss_bytes,
            tree_vms_bytes = tree$vms_bytes,
            process_count = tree$process_count,
            pids = tree$pids,
            disk_bytes = disk$total_bytes,
            file_count = disk$file_count,
            disk_category_bytes = disk$bytes
        )
        samples[[index]] <- sample
        if (file.exists(file.path(run_dir, "retention-ready")) && tree$rss_bytes > 0) {
            retained <- sample
        }
        if (!process$is_alive()) {
            break
        }
        Sys.sleep(as.numeric(execution$sampling_interval_ms) / 1000)
    }
    process$wait(timeout = 1000)
    final_disk <- baselineDirMetrics(run_dir, categories)
    rss_values <- vapply(samples, `[[`, numeric(1), "tree_rss_bytes")
    vms_values <- vapply(samples, `[[`, numeric(1), "tree_vms_bytes")
    disk_values <- vapply(samples, `[[`, numeric(1), "disk_bytes")
    process_values <- vapply(samples, `[[`, integer(1), "process_count")
    artifact_categories <- intersect(
        c("committed", "staging_snapshot", "duckdb_spill", "final"),
        names(peak_category_bytes)
    )
    artifact_disk_values <- vapply(samples, \(sample) {
        sum(unlist(sample$disk_category_bytes[artifact_categories], use.names = FALSE))
    }, numeric(1))
    if (is.null(retained)) {
        retained <- samples[[length(samples)]]
    }

    list(
        status = if (identical(process$get_exit_status(), 0L)) "passed" else "failed",
        exit_status = process$get_exit_status(),
        started_at = started_at,
        finished_at = baselineUtcNow(),
        command = basename(command),
        arguments = as.list(baselineDisplayCommandArgs(command_args, run_dir)),
        stdout = if (capture_output) file.path("final", "stdout.log") else NULL,
        stderr = if (capture_output) file.path("final", "stderr.log") else NULL,
        output_retained = capture_output,
        metrics = list(
            elapsed_seconds = proc.time()[["elapsed"]] - started,
            peak_tree_rss_bytes = max(rss_values, na.rm = TRUE),
            retained_tree_rss_bytes = retained$tree_rss_bytes,
            peak_tree_vms_bytes = max(vms_values, na.rm = TRUE),
            maximum_process_count = max(process_values, na.rm = TRUE),
            peak_disk_bytes = max(disk_values, na.rm = TRUE),
            peak_artifact_disk_bytes = max(artifact_disk_values, na.rm = TRUE),
            peak_disk_category_bytes = as.list(peak_category_bytes),
            final_disk_bytes = final_disk$total_bytes,
            final_disk_category_bytes = final_disk$bytes,
            final_file_count = final_disk$file_count,
            sample_count = length(samples)
        ),
        samples = samples
    )
}

baselineProcSelfRss <- function() {
    status_path <- sprintf("/proc/%d/status", Sys.getpid())
    if (!file.exists(status_path)) {
        return(NA_real_)
    }
    line <- grep("^VmRSS:", readLines(status_path, warn = FALSE), value = TRUE)
    if (!length(line)) {
        return(NA_real_)
    }
    as.numeric(gsub("[^0-9]", "", line[[1L]])) * 1024
}

baselineNativeMetrics <- function() {
    arrow_available <- requireNamespace("arrow", quietly = TRUE)
    duckdb_available <- requireNamespace("duckdb", quietly = TRUE)
    arrow_bytes <- if (arrow_available) {
        tryCatch(as.numeric(arrow::default_memory_pool()$bytes_allocated), error = \(...) NA_real_)
    } else {
        NA_real_
    }
    list(
        arrow = list(
            available = arrow_available,
            loaded = "arrow" %in% loadedNamespaces(),
            bytes_allocated = arrow_bytes
        ),
        duckdb = list(
            available = duckdb_available,
            loaded = "duckdb" %in% loadedNamespaces(),
            connections_opened = 0L
        )
    )
}

baselineAllocationMetrics <- function(path) {
    if (!file.exists(path)) {
        return(list(available = FALSE, allocated_bytes = NA_real_, allocation_records = 0L))
    }
    lines <- readLines(path, warn = FALSE)
    bytes <- suppressWarnings(as.numeric(sub("^([0-9]+).*$", "\\1", lines)))
    list(
        available = TRUE,
        allocated_bytes = sum(bytes[is.finite(bytes)]),
        allocation_records = sum(is.finite(bytes)),
        profile_sha256 = omicsParitySha256File(path)
    )
}

baselinePrepareFixture <- function(scenario, repo_root, run_dir) {
    if (identical(scenario$kind, "committed_fixture")) {
        fixture_path <- file.path(repo_root, scenario$fixture_path)
        if (!file.exists(fixture_path)) {
            stop(sprintf("Fixture not found: %s", scenario$fixture_path), call. = FALSE)
        }
        return(scenario$fixture_path)
    }
    if (identical(scenario$kind, "generated_scaling")) {
        source_path <- normalizePath(
            file.path(repo_root, scenario$base_fixture_path),
            mustWork = TRUE
        )
        data <- utils::read.delim(source_path, check.names = FALSE, stringsAsFactors = FALSE)
        scaled <- data[
            rep(seq_len(nrow(data)), as.integer(scenario$row_multiplier)), ,
            drop = FALSE
        ]
        fixture_path <- file.path(run_dir, "staging", "generated-scaling.tsv")
        dir.create(dirname(fixture_path), recursive = TRUE, showWarnings = FALSE)
        utils::write.table(
            scaled,
            fixture_path,
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE,
            na = ""
        )
        return(fixture_path)
    }
    if (identical(scenario$kind, "optional_private")) {
        fixture_path <- Sys.getenv(scenario$fixture_env, unset = "")
        if (!nzchar(fixture_path)) {
            return(NULL)
        }
        if (!nzchar(Sys.getenv("MULTISCHOLAR_BASELINE_FINGERPRINT_SALT", unset = ""))) {
            stop(
                "Private fixture runs require MULTISCHOLAR_BASELINE_FINGERPRINT_SALT",
                call. = FALSE
            )
        }
        return(normalizePath(fixture_path, mustWork = TRUE))
    }
    stop(sprintf("Unsupported fixture kind: %s", scenario$kind), call. = FALSE)
}

baselineSanitizePrivateSummary <- function(summary) {
    salt <- Sys.getenv("MULTISCHOLAR_BASELINE_FINGERPRINT_SALT")
    summary$input_sha256 <- digest::digest(
        paste0(salt, ":", summary$input_sha256),
        algo = "sha256",
        serialize = FALSE
    )
    summary$output_sha256 <- digest::digest(
        paste0(salt, ":", summary$output_sha256),
        algo = "sha256",
        serialize = FALSE
    )
    summary$run_count <- length(summary$runs)
    summary$protein_count <- length(summary$proteins)
    summary$peptide_count <- length(summary$peptides)
    summary$runs <- NULL
    summary$proteins <- NULL
    summary$peptides <- NULL
    summary
}

baselineBoundedQuery <- function(data, query, trace_copies = FALSE) {
    iterations <- as.integer(query$iterations)
    trace_events <- character()
    runQuery <- function() {
        iteration_timings <- numeric(iterations)
        last_selected <- NULL
        for (index in seq_len(iterations)) {
            started <- as.numeric(Sys.time())
            last_selected <- data[
                data$Run == query$run,
                unlist(query$columns, use.names = FALSE),
                drop = FALSE
            ]
            ordering <- lapply(
                last_selected[unlist(query$ordering, use.names = FALSE)],
                identity
            )
            last_selected <- last_selected[
                do.call(order, c(ordering, list(method = "radix"))), ,
                drop = FALSE
            ]
            iteration_timings[[index]] <- as.numeric(Sys.time()) - started
        }
        list(timings = iteration_timings, selected = last_selected)
    }
    if (isTRUE(trace_copies)) {
        trace_events <- capture.output(
            {
                tracemem(data)
                on.exit(untracemem(data), add = TRUE)
                query_result <- runQuery()
            },
            type = "message"
        )
    } else {
        query_result <- runQuery()
    }
    timings <- query_result$timings
    selected <- query_result$selected
    list(
        query_id = query$query_id,
        rows = nrow(selected),
        columns = ncol(selected),
        median_seconds = stats::median(timings),
        p95_seconds = omicsParityQuantile(timings, 0.95),
        maximum_seconds = max(timings),
        iterations = iterations,
        timer = "Sys.time wall-clock seconds",
        tracemem_events = sum(grepl("tracemem", trace_events, fixed = TRUE))
    )
}

baselineWorker <- function(spec_path) {
    spec <- jsonlite::read_json(spec_path, simplifyVector = FALSE)
    execution <- spec$execution
    do.call(RNGkind, as.list(unlist(execution$rng_kind, use.names = FALSE)))
    set.seed(as.integer(execution$rng_seed))
    Sys.setenv(
        TZ = execution$timezone,
        OMP_NUM_THREADS = as.character(execution$threads),
        OPENBLAS_NUM_THREADS = as.character(execution$threads),
        MKL_NUM_THREADS = as.character(execution$threads),
        ARROW_NUM_THREADS = as.character(execution$threads),
        DUCKDB_THREADS = as.character(execution$threads)
    )
    try(Sys.setlocale("LC_ALL", execution$locale), silent = TRUE)

    diagnostic_run <- isTRUE(spec$diagnostics_only)
    allocation_path <- file.path(spec$run_dir, "diagnostics", "rprofmem.out")
    if (diagnostic_run) {
        dir.create(dirname(allocation_path), recursive = TRUE, showWarnings = FALSE)
        utils::Rprofmem(allocation_path)
        on.exit(utils::Rprofmem(NULL), add = TRUE)
    }
    pkgload::load_all(
        path = spec$repo_root,
        export_all = TRUE,
        helpers = FALSE,
        attach_testthat = FALSE,
        quiet = TRUE
    )

    fixture_path <- baselinePrepareFixture(spec$scenario, spec$repo_root, spec$run_dir)
    if (is.null(fixture_path)) {
        result <- list(
            schema_version = "1.0.0",
            status = "skipped",
            reason = "private_fixture_not_configured"
        )
        baselineWriteJson(result, spec$result_path)
        return(invisible(result))
    }
    stage_started <- proc.time()[["elapsed"]]
    import_result <- suppressMessages(importDIANNData(
        fixture_path,
        use_precursor_norm = isTRUE(spec$scenario$parameters$use_precursor_norm)
    ))
    import_seconds <- proc.time()[["elapsed"]] - stage_started
    summary <- omicsParitySummarizeDiann(import_result, fixture_path, spec$comparison_contract)
    if (identical(spec$scenario$kind, "optional_private")) {
        summary <- baselineSanitizePrivateSummary(summary)
    }
    baselineWriteJson(summary, file.path(spec$run_dir, "snapshots", "import-summary.json"))

    query_results <- lapply(spec$bounded_queries, \(query) {
        baselineBoundedQuery(
            as.data.frame(import_result$data),
            query,
            trace_copies = diagnostic_run
        )
    })
    baselineWriteJson(
        query_results,
        file.path(spec$run_dir, "snapshots", "bounded-query-summary.json")
    )
    if (diagnostic_run) {
        utils::Rprofmem(NULL)
        allocation <- baselineAllocationMetrics(allocation_path)
    } else {
        allocation <- list(
            available = FALSE,
            reason = "Rprofmem and tracemem run in a separate diagnostic process"
        )
    }
    result <- list(
        schema_version = "1.0.0",
        status = "passed",
        measurement_class = if (diagnostic_run) {
            "allocation_copy_diagnostic"
        } else {
            "release_performance"
        },
        scenario_id = spec$scenario$scenario_id,
        backend = "memory",
        fixture = list(
            kind = spec$scenario$kind,
            committed_bytes = if (identical(spec$scenario$kind, "committed_fixture")) {
                as.numeric(file.info(fixture_path)$size)
            } else {
                0
            },
            fingerprint = summary$input_sha256,
            source_path_retained = FALSE
        ),
        stages = list(
            list(
                stage_id = "import",
                elapsed_seconds = import_seconds,
                retained_rss_bytes = baselineProcSelfRss()
            ),
            list(stage_id = "bounded_query", results = query_results)
        ),
        observed_summary = summary,
        allocation_diagnostics = allocation,
        thread_environment = baselineThreadEnvironment(),
        native_resources = baselineNativeMetrics(),
        retention_point = paste(
            "after import summary and bounded queries while imported table",
            "remains live"
        )
    )
    baselineWriteJson(result, spec$result_path)
    writeLines("ready", file.path(spec$run_dir, "retention-ready"))
    Sys.sleep(as.numeric(execution$retention_window_ms) / 1000)
    invisible(result)
}

baselineEnvironment <- function(manifest) {
    packages <- c("MultiScholaR", "arrow", "duckdb", "processx", "ps", "testthat")
    versions <- lapply(packages, \(package) {
        if (requireNamespace(package, quietly = TRUE)) {
            as.character(utils::packageVersion(package))
        } else {
            NULL
        }
    })
    names(versions) <- packages
    if (is.null(versions$MultiScholaR)) {
        versions$MultiScholaR <- read.dcf(
            file.path(.BASELINE_REPO_ROOT, "DESCRIPTION"),
            fields = "Version"
        )[[1L]]
    }
    git_sha <- tryCatch(
        system2(
            "git",
            c("-C", shQuote(.BASELINE_REPO_ROOT), "rev-parse", "HEAD"),
            stdout = TRUE,
            stderr = FALSE
        )[[1L]],
        error = \(...) NA_character_
    )
    memory_line <- if (file.exists("/proc/meminfo")) {
        grep("^MemTotal:", readLines("/proc/meminfo", warn = FALSE), value = TRUE)
    } else {
        character()
    }
    os_info <- Sys.info()
    public_os_fields <- intersect(c("sysname", "release", "version", "machine"), names(os_info))
    list(
        captured_at = baselineUtcNow(),
        git_sha = git_sha,
        r_version = R.version.string,
        platform = R.version$platform,
        os = as.list(os_info[public_os_fields]),
        logical_cores = parallel::detectCores(logical = TRUE),
        physical_cores = parallel::detectCores(logical = FALSE),
        total_memory = if (length(memory_line)) memory_line[[1L]] else NULL,
        locale = Sys.getlocale(),
        timezone = Sys.timezone(),
        rng_kind = as.list(RNGkind()),
        rng_seed = manifest$execution$rng_seed,
        configured_locale = manifest$execution$locale,
        configured_timezone = manifest$execution$timezone,
        configured_threads = manifest$execution$threads,
        thread_control_contract = manifest$thread_controls,
        package_versions = versions,
        parent_thread_environment = baselineThreadEnvironment()
    )
}

baselineFixtureRun <- function(
    scenario,
    repetition,
    manifest,
    output_dir,
    diagnostics_only = FALSE
) {
    run_id <- if (diagnostics_only) {
        sprintf("%s-diagnostic", scenario$scenario_id)
    } else {
        sprintf("%s-r%02d", scenario$scenario_id, repetition)
    }
    run_dir <- normalizePath(file.path(output_dir, "runs", run_id), mustWork = FALSE)
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
    result_path <- file.path(run_dir, "final", "worker-result.json")
    spec_path <- tempfile("omics-baseline-worker-", fileext = ".json")
    on.exit(unlink(spec_path), add = TRUE)
    spec <- list(
        schema_version = "1.0.0",
        repo_root = .BASELINE_REPO_ROOT,
        run_dir = run_dir,
        result_path = result_path,
        scenario = scenario,
        execution = manifest$execution,
        thread_controls = manifest$thread_controls,
        comparison_contract = manifest$comparison_contract,
        bounded_queries = manifest$bounded_queries,
        diagnostics_only = diagnostics_only
    )
    baselineWriteJson(spec, spec_path)
    measured <- baselineMeasureSubprocess(
        file.path(R.home("bin"), "Rscript"),
        c("--vanilla", baselineScriptPath(), "--worker-spec", spec_path),
        run_dir,
        manifest$execution,
        manifest$disk_categories,
        env = if (identical(scenario$kind, "optional_private")) {
            c(
                MULTISCHOLAR_DIA_BASELINE_FIXTURE = Sys.getenv(scenario$fixture_env),
                MULTISCHOLAR_BASELINE_FINGERPRINT_SALT = Sys.getenv(
                    "MULTISCHOLAR_BASELINE_FINGERPRINT_SALT"
                )
            )
        } else {
            character()
        },
        capture_output = !identical(scenario$kind, "optional_private")
    )
    worker <- if (file.exists(result_path)) {
        jsonlite::read_json(result_path, simplifyVector = FALSE)
    } else {
        list(status = "missing", reason = "worker result was not written")
    }
    measured$run_id <- run_id
    measured$scenario_id <- scenario$scenario_id
    measured$repetition <- repetition
    measured$measurement_class <- if (diagnostics_only) {
        "allocation_copy_diagnostic"
    } else {
        "release_performance"
    }
    measured$cache_state <- manifest$execution$cache_sequence[[
        min(repetition, length(manifest$execution$cache_sequence))
    ]]
    measured$worker <- worker
    measured$metrics$committed_input_bytes <- if (!is.null(worker$fixture$committed_bytes)) {
        as.numeric(worker$fixture$committed_bytes)
    } else {
        0
    }
    query_p95 <- if (!is.null(worker$stages) && length(worker$stages) >= 2L) {
        vapply(worker$stages[[2L]]$results, \(query) as.numeric(query$p95_seconds), numeric(1))
    } else {
        numeric()
    }
    measured$metrics$bounded_query_p95_seconds <- if (length(query_p95)) {
        max(query_p95)
    } else {
        NA_real_
    }

    if (!is.null(scenario$oracle_id) && identical(worker$status, "passed")) {
        oracle <- omicsParityReadOracle(file.path(.BASELINE_REPO_ROOT, manifest$oracle_path))
        matches <- Filter(
            \(entry) identical(entry$scenario_id, scenario$oracle_id),
            oracle$scenarios
        )
        if (length(matches) != 1L) {
            measured$status <- "failed"
            measured$worker$oracle_comparison <- list(
                equal = FALSE,
                errors = list("exactly one oracle entry is required")
            )
        } else {
            comparison <- omicsParityCompareSummary(
                worker$observed_summary,
                matches[[1L]]$expected,
                manifest$comparison_contract
            )
            measured$worker$oracle_comparison <- comparison
            if (!isTRUE(comparison$equal)) {
                measured$status <- "failed"
            }
        }
    }
    measured
}

baselineDeterminismChecks <- function(runs) {
    scenario_ids <- unique(vapply(runs, `[[`, character(1), "scenario_id"))
    lapply(scenario_ids, \(scenario_id) {
        scenario_runs <- Filter(\(run) identical(run$scenario_id, scenario_id), runs)
        fingerprints <- vapply(scenario_runs, \(run) {
            fingerprint <- run$worker$observed_summary$output_sha256
            if (is.null(fingerprint)) NA_character_ else fingerprint
        }, character(1))
        fingerprints <- fingerprints[!is.na(fingerprints)]
        list(
            scenario_id = scenario_id,
            repetitions = length(scenario_runs),
            observed_output_sha256 = as.list(unique(fingerprints)),
            deterministic = length(fingerprints) == length(scenario_runs) &&
                length(unique(fingerprints)) == 1L
        )
    })
}

baselineScientificRun <- function(target, manifest, output_dir) {
    run_dir <- normalizePath(
        file.path(output_dir, "scientific", target$target_id),
        mustWork = FALSE
    )
    if (identical(target$runner, "e2e_lane")) {
        command_args <- c(
            "--vanilla", "tools/ci/run-e2e-ci.R", "--lane", target$lane,
            "--browser-required", "true", "--reporter", "summary",
            "--artifact-dir", file.path(run_dir, "final", "e2e-artifacts")
        )
    } else {
        expression <- sprintf(
            "devtools::test(filter = %s, reporter = 'summary')",
            deparse(target$test_filter)
        )
        command_args <- c("--vanilla", "-e", expression)
    }
    measured <- baselineMeasureSubprocess(
        file.path(R.home("bin"), "Rscript"),
        command_args,
        run_dir,
        manifest$execution,
        manifest$disk_categories,
        env = c(
            NOT_CRAN = "true",
            MULTISCHOLAR_E2E_BROWSER_REQUIRED = "true",
            OMP_NUM_THREADS = as.character(manifest$execution$threads),
            OPENBLAS_NUM_THREADS = as.character(manifest$execution$threads),
            MKL_NUM_THREADS = as.character(manifest$execution$threads),
            ARROW_NUM_THREADS = as.character(manifest$execution$threads),
            DUCKDB_THREADS = as.character(manifest$execution$threads),
            TZ = manifest$execution$timezone
        )
    )
    measured$run_id <- target$target_id
    measured$target <- target
    measured
}

baselineMain <- function(argv = commandArgs(trailingOnly = TRUE)) {
    args <- baselineParseArgs(argv)
    if (isTRUE(args$help)) {
        baselineUsage()
        return(invisible(0L))
    }
    if (!is.null(args$worker_spec)) {
        baselineWorker(args$worker_spec)
        return(invisible(0L))
    }
    if (!args$mode %in% c("fixture", "scientific", "all")) {
        stop("--mode must be fixture, scientific, or all", call. = FALSE)
    }

    manifest_path <- baselineResolvePath(args$manifest, must_work = TRUE)
    manifest <- omicsParityReadManifest(manifest_path)
    output_dir <- baselineResolvePath(args$output_dir)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    repetitions <- if (is.null(args$repetitions)) {
        as.integer(manifest$execution$repetitions)
    } else {
        args$repetitions
    }
    if (!is.finite(repetitions) || repetitions < 1L) {
        stop("--repetitions must be a positive integer", call. = FALSE)
    }

    scenarios <- manifest$scenarios
    if (length(args$scenario)) {
        scenarios <- Filter(\(scenario) scenario$scenario_id %in% args$scenario, scenarios)
        missing <- setdiff(args$scenario, vapply(scenarios, `[[`, character(1), "scenario_id"))
        if (length(missing)) {
            stop(sprintf("Unknown scenario(s): %s", paste(missing, collapse = ", ")), call. = FALSE)
        }
    }
    if (!isTRUE(args$include_private)) {
        scenarios <- Filter(\(scenario) !identical(scenario$kind, "optional_private"), scenarios)
    }

    run_id <- sprintf(
        "omics-memory-baseline-%s-%d",
        format(Sys.time(), "%Y%m%dT%H%M%S"),
        Sys.getpid()
    )
    fixture_runs <- list()
    diagnostic_runs <- list()
    if (args$mode %in% c("fixture", "all")) {
        for (scenario in scenarios) {
            for (repetition in seq_len(repetitions)) {
                fixture_runs[[length(fixture_runs) + 1L]] <- baselineFixtureRun(
                    scenario,
                    repetition,
                    manifest,
                    output_dir
                )
            }
            if (isTRUE(args$diagnostics)) {
                diagnostic_runs[[length(diagnostic_runs) + 1L]] <- baselineFixtureRun(
                    scenario,
                    1L,
                    manifest,
                    output_dir,
                    diagnostics_only = TRUE
                )
            }
        }
    }
    scientific_runs <- list()
    if (args$mode %in% c("scientific", "all")) {
        scientific_runs <- lapply(
            manifest$scientific_targets,
            baselineScientificRun,
            manifest = manifest,
            output_dir = output_dir
        )
    }
    determinism <- baselineDeterminismChecks(fixture_runs)
    failed_determinism <- vapply(determinism, \(check) !isTRUE(check$deterministic), logical(1))
    if (any(failed_determinism)) {
        failed_ids <- vapply(determinism[failed_determinism], `[[`, character(1), "scenario_id")
        fixture_runs <- lapply(fixture_runs, \(run) {
            if (run$scenario_id %in% failed_ids) {
                run$status <- "failed"
                run$determinism_failure <- TRUE
            }
            run
        })
    }
    all_runs <- c(fixture_runs, scientific_runs)
    result <- list(
        schema_version = "1.0.0",
        run_id = run_id,
        generated_at = baselineUtcNow(),
        corpus_version = manifest$corpus_version,
        capability_id = manifest$capability_id,
        backend = "memory",
        manifest_sha256 = omicsParitySha256File(manifest_path),
        environment = baselineEnvironment(manifest),
        thread_controls = manifest$thread_controls,
        sampling = list(
            method = "ps recursive process tree sampled by the parent process",
            interval_ms = manifest$execution$sampling_interval_ms,
            process_tree = TRUE,
            child_processes_included = TRUE,
            retention_marker = "retention-ready",
            rss_is_release_metric = TRUE,
            rprofmem_is_allocation_diagnostic_only = TRUE,
            tracemem_is_copy_diagnostic_only = TRUE,
            diagnostics_use_separate_processes = TRUE,
            diagnostic_runs_excluded_from_release_summary = TRUE,
            object_size_is_release_metric = FALSE,
            garbage_collection_is_release_proof = FALSE,
            cache_control = manifest$execution$cache_control,
            concurrency = manifest$execution$concurrency
        ),
        release_gates = manifest$release_gates,
        runs = all_runs,
        diagnostics = diagnostic_runs,
        determinism = determinism,
        scenario_summaries = lapply(
            split(fixture_runs, vapply(fixture_runs, `[[`, character(1), "scenario_id")),
            omicsParitySummarizeMeasurements
        ),
        summary = omicsParitySummarizeMeasurements(all_runs)
    )
    result_path <- file.path(output_dir, "baseline-result.json")
    baselineWriteJson(result, result_path)
    run_status <- if (result$summary$failed) "failed" else "passed"
    cat(sprintf(
        "OMICS baseline result=%s runs=%d output=%s\n",
        run_status,
        length(all_runs),
        result_path
    ))
    if (result$summary$failed) {
        return(invisible(1L))
    }
    invisible(0L)
}

.BASELINE_DIRECT <- any(
    normalizePath(
        sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)),
        mustWork = FALSE
    ) == baselineScriptPath()
)
if (isTRUE(.BASELINE_DIRECT)) {
    status <- baselineMain()
    if (is.numeric(status) && status != 0L) {
        quit(status = status)
    }
}
