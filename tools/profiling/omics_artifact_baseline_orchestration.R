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
    public_os_fields <- intersect(
        c("sysname", "release", "version", "machine"),
        names(os_info)
    )
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
    backend = "memory",
    diagnostics_only = FALSE
) {
    run_id <- if (diagnostics_only) {
        sprintf("%s-%s-diagnostic", scenario$scenario_id, backend)
    } else {
        sprintf("%s-%s-r%02d", scenario$scenario_id, backend, repetition)
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
        backend = backend,
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
    if (identical(scenario$kind, "optional_private") &&
        identical(backend, "artifact")) {
        private_project <- file.path(run_dir, "committed", "project")
        unlink(private_project, recursive = TRUE, force = TRUE)
        retained <- file.exists(private_project) || dir.exists(private_project)
        worker$private_artifact_payload_retained_in_evidence <- retained
        if (isTRUE(retained)) {
            measured$status <- "failed"
            worker$status <- "failed"
            worker$reason <- "private_artifact_cleanup_failed"
        }
    }
    measured$run_id <- run_id
    measured$scenario_id <- scenario$scenario_id
    measured$backend <- backend
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
    measured$metrics$committed_input_bytes <- if (!is.null(
        worker$fixture$committed_bytes
    )) {
        as.numeric(worker$fixture$committed_bytes)
    } else {
        0
    }
    measured$metrics$committed_input_file_count <- if (!is.null(
        worker$fixture$committed_file_count
    )) {
        as.numeric(worker$fixture$committed_file_count)
    } else {
        0
    }
    measured$metrics$committed_artifact_bytes <- if (!is.null(
        worker$artifact_storage$committed_bytes
    )) {
        as.numeric(worker$artifact_storage$committed_bytes)
    } else {
        0
    }
    measured$metrics$committed_artifact_file_count <- if (!is.null(
        worker$artifact_storage$committed_file_count
    )) {
        as.numeric(worker$artifact_storage$committed_file_count)
    } else {
        0
    }
    query_stage <- Filter(function(stage) {
        identical(stage$stage_id, "bounded_query")
    }, worker$stages)
    query_p95 <- if (length(query_stage) == 1L) {
        vapply(query_stage[[1L]]$results, function(query) {
            as.numeric(query$p95_seconds)
        }, numeric(1))
    } else {
        numeric()
    }
    measured$metrics$bounded_query_p95_seconds <- if (length(query_p95)) {
        max(query_p95)
    } else {
        NA_real_
    }

    if (!is.null(scenario$oracle_id) && identical(worker$status, "passed")) {
        oracle <- omicsParityReadOracle(file.path(
            .BASELINE_REPO_ROOT,
            manifest$oracle_path
        ))
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
            if (!isTRUE(comparison$equal)) measured$status <- "failed"
        }
    }
    measured
}

baselineDeterminismChecks <- function(runs) {
    group_ids <- unique(vapply(runs, function(run) {
        paste(run$scenario_id, run$backend, sep = "::")
    }, character(1)))
    lapply(group_ids, function(group_id) {
        scenario_runs <- Filter(function(run) {
            identical(paste(run$scenario_id, run$backend, sep = "::"), group_id)
        }, runs)
        scenario_id <- scenario_runs[[1L]]$scenario_id
        fingerprints <- vapply(scenario_runs, \(run) {
            fingerprint <- run$worker$observed_summary$output_sha256
            if (is.null(fingerprint)) NA_character_ else fingerprint
        }, character(1))
        fingerprints <- fingerprints[!is.na(fingerprints)]
        list(
            scenario_id = scenario_id,
            backend = scenario_runs[[1L]]$backend,
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
