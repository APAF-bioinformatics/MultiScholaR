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
source(
    file.path(
        .BASELINE_REPO_ROOT,
        "tools",
        "profiling",
        "omics_resource_measurement.R"
    ),
    local = FALSE
)
if (!exists("omicsParityReadManifest", mode = "function")) {
    source(
        file.path(.BASELINE_REPO_ROOT, "tests", "testthat", "helper-omics-parity.R"),
        local = FALSE
    )
}
source(
    file.path(.BASELINE_REPO_ROOT, "tools", "profiling", "omics_artifact_closeout_helpers.R"),
    local = FALSE
)
source(
    file.path(
        .BASELINE_REPO_ROOT,
        "tools",
        "profiling",
        "omics_artifact_baseline_orchestration.R"
    ),
    local = FALSE
)

baselineParseBool <- function(value) {
    tolower(as.character(value)) %in% c("1", "true", "yes", "y")
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
        backend = "memory",
        scenario = character(),
        repetitions = NULL,
        diagnostics = TRUE,
        include_private = FALSE,
        require_promotion = FALSE,
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
        "  --backend <memory|artifact|paired>",
        "  --scenario <id>            Repeatable fixture scenario selector",
        "  --repetitions <n>          Override manifest repetitions",
        "  --diagnostics <bool>        Run allocation/copy diagnostics separately",
        "  --include-private <bool>   Enable the opt-in private scenario",
        "  --require-promotion <bool> Exit nonzero unless paired evidence authorizes",
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
        } else if (key %in% c(
            "include_private", "diagnostics", "require_promotion"
        )) {
            args[[key]] <- baselineParseBool(value)
        } else {
            args[[key]] <- value
        }
        index <- index + 1L
    }
    args
}

baselinePrepareFixture <- function(scenario, repo_root, run_dir) {
    if (identical(scenario$kind, "committed_fixture")) {
        fixture_path <- file.path(repo_root, scenario$fixture_path)
        if (!file.exists(fixture_path)) {
            stop(sprintf("Fixture not found: %s", scenario$fixture_path), call. = FALSE)
        }
        return(normalizePath(fixture_path, mustWork = TRUE))
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

baselineBindingEnvironment <- function(manifest) {
    packages <- c("arrow", "duckdb", "processx", "vroom", "tibble")
    versions <- lapply(packages, \(package) {
        if (requireNamespace(package, quietly = TRUE)) {
            as.character(utils::packageVersion(package))
        } else {
            NULL
        }
    })
    names(versions) <- packages
    list(
        r_version = R.version.string,
        platform = R.version$platform,
        package_version = read.dcf(
            file.path(.BASELINE_REPO_ROOT, "DESCRIPTION"),
            fields = "Version"
        )[[1L]],
        package_versions = versions,
        locale = manifest$execution$locale,
        timezone = manifest$execution$timezone,
        threads = as.integer(manifest$execution$threads),
        rng_kind = unname(unlist(manifest$execution$rng_kind, use.names = FALSE)),
        rng_seed = as.integer(manifest$execution$rng_seed)
    )
}

baselineInputFingerprint <- function(path, scenario) {
    fingerprint <- omicsParitySha256File(path)
    if (!identical(scenario$kind, "optional_private")) return(fingerprint)
    salt <- Sys.getenv("MULTISCHOLAR_BASELINE_FINGERPRINT_SALT", unset = "")
    if (!nzchar(salt)) {
        stop("Private fixture runs require a fingerprint salt", call. = FALSE)
    }
    digest::digest(
        paste0(salt, ":", fingerprint),
        algo = "sha256",
        serialize = FALSE
    )
}

baselineEvidenceBinding <- function(path, scenario, manifest) {
    binding <- list(
        schema = "multischolar.omics_artifact_evidence_binding",
        schema_version = 1L,
        input_fingerprint = baselineInputFingerprint(path, scenario),
        code_revision = baselineCodeRevision(),
        capability = list(
            capability_id = manifest$capability_id,
            omic_type = "proteomics",
            input_format = "diann",
            data_level = "peptide",
            acquisition_mode = "dia"
        ),
        parameters = scenario$parameters,
        environment = baselineBindingEnvironment(manifest)
    )
    binding$binding_sha256 <- digest::digest(
        jsonlite::toJSON(
            binding,
            auto_unbox = TRUE,
            null = "null",
            na = "null",
            digits = 17
        ),
        algo = "sha256",
        serialize = FALSE
    )
    binding
}

baselineBoundedQuery <- function(data, query, trace_copies = FALSE) {
    iterations <- as.integer(query$iterations)
    selected_run <- baselineResolveQueryRun(data, query)
    limit <- baselineQueryLimit(query)
    trace_events <- character()
    runQuery <- function() {
        iteration_timings <- numeric(iterations)
        last_selected <- NULL
        for (index in seq_len(iterations)) {
            started <- as.numeric(Sys.time())
            last_selected <- data[
                data$Run == selected_run,
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
            last_selected <- head(last_selected, limit)
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
        output_sha256 = baselineQueryDigest(selected),
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
    actual_binding <- if (
        identical(spec$worker_kind, "scientific_parity") &&
        !is.null(spec$fixture_path)
    ) {
        baselineEvidenceBinding(
            spec$fixture_path,
            spec$scenario,
            list(
                capability_id = spec$binding$capability$capability_id,
                execution = spec$execution
            )
        )
    } else {
        NULL
    }
    if (!is.null(spec$package_library)) {
        namespace <- loadNamespace(
            "MultiScholaR",
            lib.loc = spec$package_library
        )
        list2env(
            as.list(namespace, all.names = TRUE),
            envir = .GlobalEnv
        )
    } else {
        pkgload::load_all(
            path = spec$repo_root,
            export_all = TRUE,
            helpers = FALSE,
            attach_testthat = FALSE,
            quiet = TRUE
        )
    }

    fixture_path <- spec$fixture_path
    if (is.null(fixture_path)) {
        result <- list(
            schema_version = "1.0.0",
            status = "skipped",
            reason = "private_fixture_not_configured"
        )
        baselineWriteJson(result, spec$result_path)
        return(invisible(result))
    }
    if (identical(spec$worker_kind, "scientific_parity")) {
        if (!identical(
            actual_binding$binding_sha256,
            spec$binding$binding_sha256
        )) {
            result <- list(
                schema_version = "1.0.0",
                status = "failed",
                reason = "scientific_parity_binding_mismatch",
                expected_binding = spec$binding,
                observed_binding = actual_binding
            )
            baselineWriteJson(result, spec$result_path)
            return(invisible(result))
        }
        import_result <- suppressMessages(importDIANNData(
            fixture_path,
            use_precursor_norm = isTRUE(
                spec$scenario$parameters$use_precursor_norm
            )
        ))
        summary <- omicsParitySummarizeDiann(
            import_result,
            fixture_path,
            spec$comparison_contract
        )
        query_results <- lapply(spec$bounded_queries, \(query) {
            baselineBoundedQuery(as.data.frame(import_result$data), query)
        })
        hydration_digest <- artifactExactHydrationDigest(import_result$data)
        if (identical(spec$scenario$kind, "optional_private")) {
            summary <- baselineSanitizePrivateSummary(summary)
            query_results <- baselineSanitizePrivateQueries(query_results)
            hydration_digest <- baselinePrivateFingerprint(hydration_digest)
        }
        result <- list(
            schema_version = "1.0.0",
            status = "passed",
            evidence_class = "independent_scientific_parity",
            binding = spec$binding,
            observed_summary = summary,
            hydration_digest = hydration_digest,
            bounded_queries = query_results
        )
        baselineWriteJson(
            summary,
            file.path(spec$run_dir, "snapshots", "parity-import-summary.json")
        )
        baselineWriteJson(result, spec$result_path)
        return(invisible(result))
    }
    stage_started <- proc.time()[["elapsed"]]
    artifact <- NULL
    if (identical(spec$backend, "artifact")) {
        artifact <- baselineArtifactStageImport(
            fixture_path,
            spec$run_dir,
            use_precursor_norm = isTRUE(
                spec$scenario$parameters$use_precursor_norm
            )
        )
        import_result <- artifact$import_result
    } else {
        import_result <- suppressMessages(importDIANNData(
            fixture_path,
            use_precursor_norm = isTRUE(
                spec$scenario$parameters$use_precursor_norm
            )
        ))
    }
    import_seconds <- proc.time()[["elapsed"]] - stage_started
    import_shape <- list(
        rows = nrow(import_result$data),
        columns = ncol(import_result$data)
    )
    selected_runs <- lapply(spec$bounded_queries, \(query) {
        baselineResolveQueryRun(import_result$data, query)
    })
    if (identical(spec$backend, "artifact")) {
        import_result$data <- NULL
        artifact$import_result <- NULL
        import_result <- NULL
        store <- newArtifactStore(
            artifact$context$getPaths(),
            artifact$context$getIdentity()$project_id
        )
        query_session <- newArtifactQuerySession(store)
        on.exit(query_session$close(), add = TRUE)
        query_results <- Map(
            \(query, selected_run) {
                baselineArtifactBoundedQuery(
                    artifact$context,
                    artifact$ref,
                    query,
                    selected_run,
                    trace_copies = diagnostic_run,
                    query_session = query_session
                )
            },
            spec$bounded_queries,
            selected_runs
        )
        query_session$close()
    } else {
        query_results <- lapply(spec$bounded_queries, \(query) {
            baselineBoundedQuery(
                as.data.frame(import_result$data),
                query,
                trace_copies = diagnostic_run
            )
        })
    }
    if (identical(spec$scenario$kind, "optional_private")) {
        query_results <- baselineSanitizePrivateQueries(query_results)
    }
    hydration_digest <- if (identical(spec$backend, "artifact")) {
        artifact$hydration_digest
    } else {
        NULL
    }
    if (!is.null(hydration_digest) &&
        identical(spec$scenario$kind, "optional_private")) {
        hydration_digest <- baselinePrivateFingerprint(hydration_digest)
    }
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
        backend = spec$backend,
        fixture = list(
            kind = spec$scenario$kind,
            committed_bytes = as.numeric(file.info(fixture_path)$size),
            committed_file_count = 1L,
            fingerprint = spec$binding$input_fingerprint,
            source_path_retained = FALSE
        ),
        binding = spec$binding,
        stages = list(
            list(
                stage_id = "import",
                elapsed_seconds = import_seconds,
                retained_rss_bytes = baselineProcSelfRss()
            ),
            if (identical(spec$backend, "artifact")) list(
                stage_id = "artifact_stage_worker",
                elapsed_seconds = import_seconds,
                process_evidence = artifact$process_evidence,
                generation_id_retained = FALSE
            ) else NULL,
            list(stage_id = "bounded_query", results = query_results)
        ),
        workflow_evidence = list(
            import_shape = import_shape,
            hydration_digest = hydration_digest
        ),
        allocation_diagnostics = allocation,
        thread_environment = baselineThreadEnvironment(),
        native_resources = baselineNativeMetrics(),
        retention_point = if (identical(spec$backend, "artifact")) {
            "after artifact persistence, source-table eviction, and bounded queries"
        } else {
            paste(
                "after import summary and bounded queries while imported table",
                "remains live"
            )
        }
    )
    result$stages <- Filter(Negate(is.null), result$stages)
    if (identical(spec$backend, "artifact")) {
        result$artifact_storage <- baselineArtifactProjectMetrics(artifact$project_root)
        result$resource_evidence <- baselineArtifactResourceEvidence(artifact$project_root)
        result$source_table_retained <- FALSE
        if (identical(spec$scenario$kind, "optional_private")) {
            result$private_artifact_payload_retained_at_worker_exit <- TRUE
        }
    } else {
        result$source_table_retained <- TRUE
    }
    baselineWriteJson(result, spec$result_path)
    writeLines("ready", file.path(spec$run_dir, "retention-ready"))
    Sys.sleep(as.numeric(execution$retention_window_ms) / 1000)
    invisible(result)
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
    if (!args$backend %in% c("memory", "artifact", "paired")) {
        stop("--backend must be memory, artifact, or paired", call. = FALSE)
    }
    if (!identical(args$backend, "memory") && !identical(args$mode, "fixture")) {
        stop("artifact and paired benchmarks require --mode fixture", call. = FALSE)
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
    package_library <- if (
        args$mode %in% c("fixture", "all") && length(scenarios) > 0L
    ) {
        baselineInstallMeasurementPackage(output_dir)
    } else {
        NULL
    }
    if (!is.null(package_library)) {
        on.exit(unlink(package_library, recursive = TRUE, force = TRUE), add = TRUE)
    }

    run_id <- sprintf(
        "omics-%s-baseline-%s-%d",
        args$backend,
        format(Sys.time(), "%Y%m%dT%H%M%S"),
        Sys.getpid()
    )
    fixture_runs <- list()
    diagnostic_runs <- list()
    if (args$mode %in% c("fixture", "all")) {
        for (scenario in scenarios) {
            backends <- if (identical(args$backend, "paired")) {
                c("memory", "artifact")
            } else {
                args$backend
            }
            for (repetition in seq_len(repetitions)) {
                for (backend in backends) {
                    fixture_runs[[length(fixture_runs) + 1L]] <- baselineFixtureRun(
                        scenario,
                        repetition,
                        manifest,
                        output_dir,
                        backend = backend,
                        package_library = package_library
                    )
                }
            }
            if (isTRUE(args$diagnostics)) {
                for (backend in backends) {
                    diagnostic_runs[[length(diagnostic_runs) + 1L]] <- baselineFixtureRun(
                        scenario,
                        1L,
                        manifest,
                        output_dir,
                        backend = backend,
                        diagnostics_only = TRUE,
                        package_library = package_library
                    )
                }
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
        failed_groups <- vapply(determinism[failed_determinism], function(check) {
            paste(check$scenario_id, check$backend, sep = "::")
        }, character(1))
        fixture_runs <- lapply(fixture_runs, \(run) {
            group_id <- paste(run$scenario_id, run$backend, sep = "::")
            if (group_id %in% failed_groups) {
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
        backend = args$backend,
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
            split(fixture_runs, vapply(fixture_runs, function(run) {
                if (identical(args$backend, "memory")) {
                    run$scenario_id
                } else {
                    paste(run$scenario_id, run$backend, sep = "::")
                }
            }, character(1))),
            omicsParitySummarizeMeasurements
        ),
        summary = omicsParitySummarizeMeasurements(all_runs)
    )
    if (identical(args$backend, "paired")) {
        result$promotion_evaluation <- baselinePromotionEvaluation(
            fixture_runs,
            scenarios,
            manifest$release_gates
        )
    }
    result$promotion_required <- isTRUE(args$require_promotion)
    result_path <- file.path(output_dir, "baseline-result.json")
    baselineWriteJson(result, result_path)
    promotion_failed <- isTRUE(args$require_promotion) &&
        (!identical(args$backend, "paired") ||
            !isTRUE(result$promotion_evaluation$authorized))
    run_status <- if (result$summary$failed || promotion_failed) {
        "failed"
    } else {
        "passed"
    }
    cat(sprintf(
        "OMICS baseline result=%s runs=%d output=%s\n",
        run_status,
        length(all_runs),
        result_path
    ))
    if (result$summary$failed || promotion_failed) {
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
