baselineQueryLimit <- function(query) {
    limit <- query$limit
    if (is.null(limit)) limit <- 10000L
    limit <- as.integer(limit)
    if (!is.finite(limit) || limit < 1L) {
        stop("bounded query limit must be a positive integer", call. = FALSE)
    }
    limit
}

baselineResolveQueryRun <- function(data, query) {
    runs <- sort(unique(as.character(data$Run)))
    if (!length(runs)) {
        stop("DIA-NN bounded query requires at least one run", call. = FALSE)
    }
    requested <- as.character(query$run)
    if (length(requested) == 1L && !is.na(requested) && requested %in% runs) {
        return(requested)
    }
    runs[[1L]]
}

baselineQueryDigest <- function(data) {
    keys <- intersect(
        c("Run", "Protein.Group", "Stripped.Sequence", "Precursor.Id"),
        names(data)
    )
    omicsParityTableSha256(data, keys, names(data))
}

baselinePrivateFingerprint <- function(fingerprint) {
    digest::digest(
        paste0(Sys.getenv("MULTISCHOLAR_BASELINE_FINGERPRINT_SALT"), ":", fingerprint),
        algo = "sha256",
        serialize = FALSE
    )
}

baselineSanitizePrivateQueries <- function(results) {
    lapply(results, function(result) {
        result$output_sha256 <- baselinePrivateFingerprint(result$output_sha256)
        result
    })
}

baselineArtifactProjectMetrics <- function(project_root) {
    files <- list.files(
        project_root,
        recursive = TRUE,
        full.names = TRUE,
        all.files = TRUE,
        no.. = TRUE
    )
    files <- files[file.exists(files) & !dir.exists(files)]
    list(
        committed_bytes = sum(as.numeric(file.info(files)$size), na.rm = TRUE),
        committed_file_count = length(files)
    )
}

baselineArtifactPersistImport <- function(
    import_result,
    fixture_path,
    run_dir,
    use_precursor_norm
) {
    project_root <- file.path(run_dir, "committed", "project")
    dir.create(project_root, recursive = TRUE, showWarnings = FALSE)
    paths <- list(
        base_dir = project_root,
        project_id = "dia-closeout-project",
        omic_type = "proteomics",
        omic_label = "dia-closeout",
        source_dir = dirname(fixture_path),
        results_dir = file.path(project_root, "results")
    )
    workflow <- new.env(parent = emptyenv())
    workflow$workflow_context <- createWorkflowContext(
        paths,
        "proteomics",
        "dia-closeout",
        storage_policy = list(
            requested_backend = "artifact",
            requested_rollout = "dual_write",
            migration_requested = TRUE,
            project_id = paths$project_id
        )
    )
    workflow$data_tbl <- import_result$data
    workflow$data_format <- "diann"
    workflow$data_type <- import_result$data_type
    workflow$artifact_stage_results <- NULL
    persisted <- persistProtDiaImportArtifacts(
        workflow,
        import_result,
        fixture_path,
        use_precursor_norm = isTRUE(use_precursor_norm),
        retain_source_uri = FALSE,
        log_warn = function(...) NULL
    )
    if (!isTRUE(persisted$enabled) || !isTRUE(persisted$ok) ||
        !isTRUE(persisted$committed)) {
        stop("DIA-NN artifact benchmark import did not commit", call. = FALSE)
    }
    workflow$data_tbl <- NULL
    list(
        context = workflow$workflow_context,
        ref = persisted$refs$canonical_data,
        project_root = project_root,
        generation_id = persisted$generation_id
    )
}

baselineArtifactBoundedQuery <- function(
    context,
    ref,
    query,
    selected_run,
    trace_copies = FALSE
) {
    descriptor <- artifactDiaWorkflowDescriptor()
    operation_id <- names(descriptor$queries)[[1L]]
    iterations <- as.integer(query$iterations)
    limit <- baselineQueryLimit(query)
    timings <- numeric(iterations)
    selected <- NULL
    trace_events <- character()
    runQuery <- function() {
        for (index in seq_len(iterations)) {
            started <- as.numeric(Sys.time())
            selected <- previewProtDiaImportArtifact(
                context,
                ref,
                projections = unlist(query$columns, use.names = FALSE),
                filters = list(run = list(operator = "equal", value = selected_run)),
                limit = limit
            )
            timings[[index]] <- as.numeric(Sys.time()) - started
        }
        list(timings = timings, selected = selected)
    }
    if (isTRUE(trace_copies)) {
        trace_events <- capture.output(
            query_result <- runQuery(),
            type = "message"
        )
    } else {
        query_result <- runQuery()
    }
    timings <- query_result$timings
    selected <- query_result$selected
    list(
        query_id = query$query_id,
        operation_id = operation_id,
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

baselineArtifactResourceEvidence <- function(project_root) {
    temporary <- file.path(project_root, "state", "tmp", "duckdb")
    data_root <- file.path(project_root, "data")
    staging_dirs <- if (dir.exists(data_root)) {
        list.dirs(data_root, recursive = TRUE, full.names = TRUE)
    } else {
        character()
    }
    staging_dirs <- staging_dirs[basename(staging_dirs) == "staging"]
    staging_entries <- sum(vapply(staging_dirs, function(path) {
        length(list.files(path, all.files = TRUE, no.. = TRUE))
    }, integer(1)))
    list(
        open_r_connections = nrow(showConnections(all = TRUE)),
        duckdb_temporary_entries = if (dir.exists(temporary)) {
            length(list.files(temporary, all.files = TRUE, no.. = TRUE))
        } else {
            0L
        },
        staging_entries = staging_entries
    )
}

baselineMetricValue <- function(run, field) {
    value <- run$metrics[[field]]
    if (is.null(value) || length(value) != 1L) return(NA_real_)
    as.numeric(value)
}

baselinePairedValues <- function(pairs, memory_field, artifact_field, operation) {
    values <- vapply(pairs, function(pair) {
        memory <- baselineMetricValue(pair$memory, memory_field)
        artifact <- baselineMetricValue(pair$artifact, artifact_field)
        if (!is.finite(memory) || memory <= 0 || !is.finite(artifact)) return(NA_real_)
        operation(memory, artifact)
    }, numeric(1))
    values[is.finite(values)]
}

baselineGateEvidence <- function(gate_id, values, threshold, comparison) {
    values <- as.numeric(values)
    values <- values[is.finite(values)]
    observed <- if (length(values)) stats::median(values) else NA_real_
    passed <- is.finite(observed) && switch(
        comparison,
        minimum = observed >= threshold,
        maximum = observed <= threshold,
        stop("unsupported gate comparison", call. = FALSE)
    )
    list(
        gate_id = gate_id,
        comparison = comparison,
        threshold = threshold,
        observed_median = observed,
        observed_p95 = omicsParityQuantile(values, 0.95),
        valid_pairs = length(values),
        passed = passed
    )
}

baselinePairScenarioRuns <- function(runs, scenario_id) {
    selected <- Filter(
        function(run) identical(run$scenario_id, scenario_id) &&
            identical(run$status, "passed") &&
            identical(run$worker$status, "passed"),
        runs
    )
    repetitions <- sort(unique(vapply(selected, `[[`, integer(1), "repetition")))
    pairs <- lapply(repetitions, function(repetition) {
        candidates <- Filter(
            function(run) identical(run$repetition, repetition),
            selected
        )
        memory <- Filter(function(run) identical(run$backend, "memory"), candidates)
        artifact <- Filter(function(run) identical(run$backend, "artifact"), candidates)
        if (length(memory) != 1L || length(artifact) != 1L) return(NULL)
        list(memory = memory[[1L]], artifact = artifact[[1L]])
    })
    Filter(Negate(is.null), pairs)
}

baselinePairScientificParity <- function(pairs) {
    queryResults <- function(worker) {
        stages <- Filter(function(stage) {
            identical(stage$stage_id, "bounded_query")
        }, worker$stages)
        if (length(stages) != 1L) return(NULL)
        vapply(stages[[1L]]$results, `[[`, character(1), "output_sha256")
    }
    vapply(pairs, function(pair) {
        memory <- pair$memory$worker
        artifact <- pair$artifact$worker
        summaries_equal <- identical(
            memory$observed_summary$output_sha256,
            artifact$observed_summary$output_sha256
        )
        memory_queries <- queryResults(memory)
        artifact_queries <- queryResults(artifact)
        summaries_equal && !is.null(memory_queries) &&
            identical(memory_queries, artifact_queries)
    }, logical(1))
}

baselineScenarioPromotionEvaluation <- function(scenario, runs, gates) {
    pairs <- baselinePairScenarioRuns(runs, scenario$scenario_id)
    ratio <- function(memory, artifact) artifact / memory
    reduction <- function(memory, artifact) 1 - artifact / memory
    gate_results <- list(
        baselineGateEvidence(
            "minimum_peak_rss_reduction_fraction",
            baselinePairedValues(
                pairs, "peak_tree_rss_bytes", "peak_tree_rss_bytes", reduction
            ),
            gates$minimum_peak_rss_reduction_fraction,
            "minimum"
        ),
        baselineGateEvidence(
            "minimum_retained_rss_reduction_fraction",
            baselinePairedValues(
                pairs, "retained_tree_rss_bytes", "retained_tree_rss_bytes", reduction
            ),
            gates$minimum_retained_rss_reduction_fraction,
            "minimum"
        ),
        baselineGateEvidence(
            "maximum_runtime_ratio",
            baselinePairedValues(pairs, "elapsed_seconds", "elapsed_seconds", ratio),
            gates$maximum_runtime_ratio,
            "maximum"
        ),
        baselineGateEvidence(
            "maximum_committed_disk_ratio",
            baselinePairedValues(
                pairs, "committed_input_bytes", "committed_artifact_bytes", ratio
            ),
            gates$maximum_committed_disk_ratio,
            "maximum"
        ),
        baselineGateEvidence(
            "maximum_peak_disk_ratio",
            baselinePairedValues(
                pairs, "committed_input_bytes", "peak_artifact_disk_bytes", ratio
            ),
            gates$maximum_peak_disk_ratio,
            "maximum"
        ),
        baselineGateEvidence(
            "maximum_final_file_count_ratio",
            baselinePairedValues(
                pairs,
                "final_file_count",
                "final_file_count",
                ratio
            ),
            gates$maximum_final_file_count_ratio,
            "maximum"
        ),
        baselineGateEvidence(
            "maximum_bounded_query_p95_ratio",
            baselinePairedValues(
                pairs,
                "bounded_query_p95_seconds",
                "bounded_query_p95_seconds",
                ratio
            ),
            gates$maximum_bounded_query_p95_ratio,
            "maximum"
        ),
        baselineGateEvidence(
            "maximum_bounded_query_p95_seconds",
            vapply(pairs, function(pair) {
                baselineMetricValue(pair$artifact, "bounded_query_p95_seconds")
            }, numeric(1)),
            gates$maximum_bounded_query_p95_seconds,
            "maximum"
        )
    )
    minimum_pairs <- as.integer(gates$minimum_representative_repetitions)
    gate_results <- lapply(gate_results, function(gate) {
        gate$required_pairs <- minimum_pairs
        gate$passed <- isTRUE(gate$passed) && gate$valid_pairs >= minimum_pairs
        gate
    })
    parity <- baselinePairScientificParity(pairs)
    representative <- identical(scenario$kind, "optional_private")
    enough_pairs <- length(pairs) >= minimum_pairs
    all_gates <- length(gate_results) > 0L && all(vapply(
        gate_results,
        `[[`,
        logical(1),
        "passed"
    ))
    parity_passed <- length(parity) == length(pairs) && length(parity) > 0L && all(parity)
    list(
        scenario_id = scenario$scenario_id,
        fixture_kind = scenario$kind,
        promotion_eligible = representative,
        paired_repetitions = length(pairs),
        minimum_representative_repetitions = gates$minimum_representative_repetitions,
        scientific_parity = list(
            valid_pairs = length(parity),
            passed = parity_passed
        ),
        gates = gate_results,
        authorized = representative && enough_pairs && parity_passed && all_gates,
        reason = if (!representative) {
            "non_representative_evidence"
        } else if (!enough_pairs) {
            "insufficient_paired_repetitions"
        } else if (!parity_passed) {
            "scientific_or_query_parity_failed"
        } else if (!all_gates) {
            "one_or_more_performance_gates_failed"
        } else {
            "all_representative_gates_passed"
        }
    )
}

baselinePromotionEvaluation <- function(runs, scenarios, gates) {
    evaluations <- lapply(
        scenarios,
        baselineScenarioPromotionEvaluation,
        runs = runs,
        gates = gates
    )
    candidates <- Filter(function(result) isTRUE(result$promotion_eligible), evaluations)
    authorized <- length(candidates) > 0L && all(vapply(
        candidates,
        `[[`,
        logical(1),
        "authorized"
    ))
    list(
        schema_version = "1.0.0",
        capability_id = "proteomics.diann.peptide.dia.v1",
        authorization_scope = "exact_dia_nn_peptide_dia_tuple_only",
        generated_evidence_authorizes_promotion = FALSE,
        representative_candidates = length(candidates),
        authorized = authorized,
        reason = if (!length(candidates)) {
            "representative_scenario_not_run"
        } else if (!authorized) {
            "representative_candidate_failed"
        } else {
            "all_representative_candidates_passed"
        },
        scenarios = evaluations
    )
}

baselineValidateCloseoutResult <- function(path) {
    result <- omicsParityReadJson(
        path,
        .OMICS_PARITY_RESULT_SCHEMA,
        "DIA-NN artifact closeout result"
    )
    required <- c(
        "run_id", "backend", "capability_id", "environment", "runs",
        "release_gates", "promotion_evaluation", "summary"
    )
    missing <- setdiff(required, names(result))
    if (length(missing)) {
        stop(
            sprintf("Closeout result is missing: %s", paste(missing, collapse = ", ")),
            call. = FALSE
        )
    }
    valid <- identical(result$backend, "paired") &&
        identical(result$capability_id, "proteomics.diann.peptide.dia.v1") &&
        identical(result$promotion_evaluation$schema_version, "1.0.0") &&
        identical(
            result$promotion_evaluation$authorization_scope,
            "exact_dia_nn_peptide_dia_tuple_only"
        ) &&
        is.logical(result$promotion_evaluation$authorized) &&
        length(result$promotion_evaluation$authorized) == 1L &&
        all(vapply(result$runs, function(run) {
            run$backend %in% c("memory", "artifact")
        }, logical(1)))
    if (!isTRUE(valid)) {
        stop("DIA-NN artifact closeout result is incompatible", call. = FALSE)
    }
    result
}
