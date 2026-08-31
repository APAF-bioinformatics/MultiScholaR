publicationStatisticsAbort <- function(message) {
    publicationAbort(message, "multischolar_publication_statistics_error")
}

publicationPairMetricTable <- function(
    records,
    metric,
    numerator_arm = "B",
    denominator_arm = "A"
) {
    required <- c("project_id", "pair_id", "arm", "status", metric)
    if (!is.data.frame(records) || !all(required %in% names(records)) ||
        !numerator_arm %in% c("A", "B") || !denominator_arm %in% c("A", "B") ||
        identical(numerator_arm, denominator_arm)) {
        publicationStatisticsAbort("Pair metric records are malformed")
    }
    records <- records[records$status == "passed", , drop = FALSE]
    if (!nrow(records)) {
        publicationStatisticsAbort("No passed pairs are available")
    }
    binding_fields <- intersect(
        c("source_id", "session_id", "block_id", "host_id"),
        names(records)
    )
    pair_ids <- unique(records$pair_id)
    pairs <- lapply(pair_ids, \(pair_id) {
        rows <- records[records$pair_id == pair_id, , drop = FALSE]
        numerator <- rows[rows$arm == numerator_arm, , drop = FALSE]
        denominator <- rows[rows$arm == denominator_arm, , drop = FALSE]
        if (nrow(numerator) != 1L || nrow(denominator) != 1L ||
            !identical(numerator$project_id, denominator$project_id)) {
            publicationStatisticsAbort("Pair arm or project ownership differs")
        }
        if (any(!vapply(binding_fields, \(field) {
            identical(numerator[[field]], denominator[[field]])
        }, logical(1)))) {
            publicationStatisticsAbort("Pair source session block or host differs")
        }
        numerator_value <- as.numeric(numerator[[metric]])
        denominator_value <- as.numeric(denominator[[metric]])
        if (!publicationScalarNumber(numerator_value, positive = TRUE) ||
            !publicationScalarNumber(denominator_value, positive = TRUE)) {
            publicationStatisticsAbort("Pair metric must be finite and positive")
        }
        row <- data.frame(
            project_id = numerator$project_id,
            pair_id = pair_id,
            numerator = numerator_value,
            denominator = denominator_value,
            ratio = numerator_value / denominator_value,
            log_ratio = log(numerator_value / denominator_value),
            difference = numerator_value - denominator_value,
            stringsAsFactors = FALSE
        )
        for (field in binding_fields) row[[field]] <- numerator[[field]]
        row
    })
    result <- do.call(rbind, pairs)
    rownames(result) <- NULL
    result
}

publicationQuantile <- function(values, probability) {
    as.numeric(stats::quantile(
        values,
        probs = probability,
        names = FALSE,
        type = 8,
        na.rm = FALSE
    ))
}

publicationRatioSummary <- function(pairs) {
    if (!is.data.frame(pairs) || !nrow(pairs) ||
        any(!is.finite(pairs$log_ratio))) {
        publicationStatisticsAbort("Pair table cannot be summarized")
    }
    list(
        pairs = nrow(pairs),
        geometric_mean_ratio = exp(mean(pairs$log_ratio)),
        median_ratio = stats::median(pairs$ratio),
        median_difference = stats::median(pairs$difference),
        median_percent_change = 100 * (stats::median(pairs$ratio) - 1),
        ratio_iqr = unname(as.list(stats::quantile(
            pairs$ratio,
            c(0.25, 0.75),
            names = FALSE,
            type = 8
        ))),
        ratio_minimum = min(pairs$ratio),
        ratio_maximum = max(pairs$ratio)
    )
}

publicationBootstrapWithinProject <- function(
    pairs,
    resamples = 10000L,
    seed = 106231L,
    confidence = 0.95
) {
    if (!is.data.frame(pairs) || !nrow(pairs) || resamples < 1000L) {
        publicationStatisticsAbort("Within-project bootstrap inputs are invalid")
    }
    estimates <- publicationWithPreservedSeed(seed, replicate(resamples, {
        sampled <- pairs[sample.int(nrow(pairs), replace = TRUE), , drop = FALSE]
        mean(sampled$log_ratio)
    }))
    alpha <- (1 - confidence) / 2
    list(
        resampling_unit = "process_pair",
        resamples = as.integer(resamples),
        seed = as.integer(seed),
        confidence = confidence,
        geometric_mean_ratio = exp(mean(pairs$log_ratio)),
        lower = exp(publicationQuantile(estimates, alpha)),
        upper = exp(publicationQuantile(estimates, 1 - alpha))
    )
}

publicationBootstrapHierarchical <- function(
    pairs,
    resamples = 10000L,
    seed = 106231L,
    confidence = 0.95
) {
    projects <- unique(pairs$project_id)
    if (!is.data.frame(pairs) || length(projects) < 3L || resamples < 1000L) {
        publicationStatisticsAbort("Hierarchical bootstrap requires three projects")
    }
    by_project <- split(pairs, pairs$project_id)
    estimates <- publicationWithPreservedSeed(seed, replicate(resamples, {
        sampled_projects <- sample(projects, length(projects), replace = TRUE)
        sampled_project_effects <- vapply(sampled_projects, \(project_id) {
            project_pairs <- by_project[[project_id]]
            sampled <- project_pairs[
                sample.int(nrow(project_pairs), replace = TRUE),
                ,
                drop = FALSE
            ]
            mean(sampled$log_ratio)
        }, numeric(1))
        mean(sampled_project_effects)
    }))
    alpha <- (1 - confidence) / 2
    project_effects <- lapply(by_project, publicationRatioSummary)
    project_log_means <- vapply(
        by_project,
        \(project_pairs) mean(project_pairs$log_ratio),
        numeric(1)
    )
    estimate <- mean(project_log_means)
    list(
        resampling_unit = "project_then_process_pair",
        projects = length(projects),
        pairs = nrow(pairs),
        resamples = as.integer(resamples),
        seed = as.integer(seed),
        confidence = confidence,
        project_weighting = "equal_project_weight_after_within_project_resampling",
        geometric_mean_ratio = exp(estimate),
        lower = exp(publicationQuantile(estimates, alpha)),
        upper = exp(publicationQuantile(estimates, 1 - alpha)),
        project_effects = project_effects
    )
}

publicationSummarizeMetric <- function(
    records,
    metric,
    numerator_arm = "B",
    denominator_arm = "A",
    resamples = 10000L,
    seed = 106231L
) {
    source_result_digest <- publicationObjectDigest(as.list(records))
    passed <- records[records$status == "passed", , drop = FALSE]
    pair_ids <- unique(as.character(records$pair_id))
    complete_pair_ids <- pair_ids[vapply(pair_ids, \(pair_id) {
        pair <- passed[passed$pair_id == pair_id, , drop = FALSE]
        sum(pair$arm == numerator_arm) == 1L &&
            sum(pair$arm == denominator_arm) == 1L
    }, logical(1))]
    excluded_pair_ids <- setdiff(pair_ids, complete_pair_ids)
    malformed_pair_ids <- pair_ids[vapply(pair_ids, \(pair_id) {
        pair <- records[records$pair_id == pair_id, , drop = FALSE]
        passed_pair <- pair[pair$status == "passed", , drop = FALSE]
        sum(passed_pair$arm == numerator_arm) > 1L ||
            sum(passed_pair$arm == denominator_arm) > 1L ||
            (!any(pair$status != "passed") && !pair_id %in% complete_pair_ids)
    }, logical(1))]
    if (length(malformed_pair_ids)) {
        publicationStatisticsAbort("Malformed pair records cannot be excluded")
    }
    analysis_records <- passed[
        passed$pair_id %in% complete_pair_ids,
        ,
        drop = FALSE
    ]
    pairs <- publicationPairMetricTable(
        analysis_records,
        metric,
        numerator_arm,
        denominator_arm
    )
    by_project <- lapply(split(pairs, pairs$project_id), \(project_pairs) {
        list(
            summary = publicationRatioSummary(project_pairs),
            bootstrap = publicationBootstrapWithinProject(
                project_pairs,
                resamples,
                seed
            )
        )
    })
    session_results <- if ("session_id" %in% names(pairs)) {
        lapply(split(pairs, pairs$session_id), publicationRatioSummary)
    } else {
        list()
    }
    list(
        metric = metric,
        numerator_arm = numerator_arm,
        denominator_arm = denominator_arm,
        analysis_seed = as.integer(seed),
        bootstrap_resamples = as.integer(resamples),
        source_result_digest = source_result_digest,
        pair_digest = publicationObjectDigest(as.list(pairs)),
        all_observations = records,
        excluded_pair_ids = as.list(excluded_pair_ids),
        excluded_pair_count = length(excluded_pair_ids),
        paired_observations = pairs,
        failure_summary = publicationFailureSummary(records),
        arm_variability = list(
            numerator = publicationArmVariability(
                records,
                metric,
                numerator_arm
            ),
            denominator = publicationArmVariability(
                records,
                metric,
                denominator_arm
            )
        ),
        summary = publicationRatioSummary(pairs),
        project_results = by_project,
        session_results = session_results,
        hierarchical = publicationBootstrapHierarchical(pairs, resamples, seed)
    )
}

publicationArmVariability <- function(records, metric, arm) {
    values <- as.numeric(records[
        records$arm == arm & records$status == "passed",
        metric,
        drop = TRUE
    ])
    if (!length(values) || any(!is.finite(values)) || any(values <= 0)) {
        publicationStatisticsAbort("Arm variability values are invalid")
    }
    quartiles <- stats::quantile(values, c(0.25, 0.75), names = FALSE, type = 8)
    list(
        arm = arm,
        observations = length(values),
        median = stats::median(values),
        q25 = quartiles[[1L]],
        q75 = quartiles[[2L]],
        minimum = min(values),
        maximum = max(values),
        mean = mean(values),
        standard_deviation = if (length(values) > 1L) stats::sd(values) else 0,
        coefficient_of_variation = if (length(values) > 1L) {
            stats::sd(values) / mean(values)
        } else {
            0
        }
    )
}

publicationFailureSummary <- function(records) {
    statuses <- sort(unique(as.character(records$status)), method = "radix")
    counts <- stats::setNames(lapply(statuses, \(status) {
        sum(records$status == status)
    }), statuses)
    failed_pair_ids <- if ("pair_id" %in% names(records)) {
        unique(as.character(records$pair_id[records$status != "passed"]))
    } else {
        character()
    }
    list(
        observations = nrow(records),
        passed = sum(records$status == "passed"),
        failed = sum(records$status != "passed"),
        failed_pairs = length(failed_pair_ids),
        failed_pair_ids = as.list(failed_pair_ids),
        timed_out = if ("timed_out" %in% names(records)) {
            sum(records$timed_out %in% TRUE, na.rm = TRUE)
        } else {
            0L
        },
        by_status = counts
    )
}

publicationQueryRunTable <- function(query_records, minimum_queries = 200L) {
    required <- c(
        "project_id", "pair_id", "arm", "run_id", "status", "query_id",
        "latency_seconds", "query_suite_elapsed_seconds", "returned_rows",
        "returned_bytes", "output_digest"
    )
    if (!is.data.frame(query_records) || !nrow(query_records) ||
        !all(required %in% names(query_records)) ||
        !publicationScalarNumber(minimum_queries, positive = TRUE)) {
        publicationStatisticsAbort("Query records are malformed")
    }
    passed <- query_records$status == "passed"
    digest_valid <- !is.na(query_records$output_digest) &
        grepl("^[0-9a-f]{64}$", query_records$output_digest)
    numeric_valid <- vapply(seq_len(nrow(query_records)), \(index) {
        if (!passed[[index]]) return(TRUE)
        publicationScalarNumber(
            query_records$latency_seconds[[index]],
            positive = TRUE
        ) && publicationScalarNumber(query_records$returned_rows[[index]]) &&
            query_records$returned_rows[[index]] >= 0 &&
            publicationScalarNumber(query_records$returned_bytes[[index]]) &&
            query_records$returned_bytes[[index]] >= 0
    }, logical(1))
    if (any(passed & (!digest_valid | !numeric_valid))) {
        publicationStatisticsAbort("Passed query evidence is invalid")
    }
    binding_fields <- intersect(
        c("project_id", "pair_id", "arm", "source_id", "session_id",
            "block_id", "host_id"),
        names(query_records)
    )
    by_run <- split(query_records, query_records$run_id)
    rows <- lapply(names(by_run), \(run_id) {
        run <- by_run[[run_id]]
        if (any(vapply(binding_fields, \(field) {
            length(unique(run[[field]])) != 1L
        }, logical(1))) || anyDuplicated(run$query_id)) {
            publicationStatisticsAbort("Query run bindings or ids differ")
        }
        run_passed <- all(run$status == "passed") && nrow(run) >= minimum_queries
        latencies <- if (run_passed) run$latency_seconds else NA_real_
        suite_elapsed <- unique(run$query_suite_elapsed_seconds)
        suite_valid <- length(suite_elapsed) == 1L &&
            publicationScalarNumber(suite_elapsed, positive = TRUE) &&
            (!run_passed || suite_elapsed >= sum(latencies))
        if (!suite_valid) {
            publicationStatisticsAbort("Query suite elapsed boundary is invalid")
        }
        row <- data.frame(
            project_id = run$project_id[[1L]],
            pair_id = run$pair_id[[1L]],
            arm = run$arm[[1L]],
            run_id = run_id,
            status = if (run_passed) "passed" else "query_contract_failed",
            query_count = nrow(run),
            query_p50_seconds = if (run_passed) {
                publicationQuantile(latencies, 0.50)
            } else {
                NA_real_
            },
            query_p95_seconds = if (run_passed) {
                publicationQuantile(latencies, 0.95)
            } else {
                NA_real_
            },
            query_p99_seconds = if (run_passed) {
                publicationQuantile(latencies, 0.99)
            } else {
                NA_real_
            },
            query_throughput_per_second = if (run_passed) {
                nrow(run) / suite_elapsed
            } else {
                NA_real_
            },
            returned_rows = sum(
                run$returned_rows[run$status == "passed"],
                na.rm = TRUE
            ),
            returned_bytes = sum(
                run$returned_bytes[run$status == "passed"],
                na.rm = TRUE
            ),
            query_result_digest = publicationObjectDigest(as.list(run)),
            stringsAsFactors = FALSE
        )
        for (field in setdiff(binding_fields, c("project_id", "pair_id", "arm"))) {
            row[[field]] <- run[[field]][[1L]]
        }
        row
    })
    result <- do.call(rbind, rows)
    rownames(result) <- NULL
    result
}

publicationSummarizeClusteredQueries <- function(
    query_records,
    minimum_queries = 200L,
    resamples = 10000L,
    seed = 106231L
) {
    runs <- publicationQueryRunTable(query_records, minimum_queries)
    metrics <- c(
        "query_p50_seconds", "query_p95_seconds", "query_p99_seconds",
        "query_throughput_per_second"
    )
    list(
        resampling_unit = "fresh_process_run_pair",
        within_run_queries_are_independent_replicates = FALSE,
        minimum_measured_queries_per_run = as.integer(minimum_queries),
        analysis_seed = as.integer(seed),
        source_result_digest = publicationObjectDigest(as.list(query_records)),
        all_query_observations = query_records,
        run_observations = runs,
        metric_results = stats::setNames(lapply(seq_along(metrics), \(index) {
            publicationSummarizeMetric(
                runs,
                metrics[[index]],
                resamples = resamples,
                seed = seed + index - 1L
            )
        }), metrics)
    )
}

publicationHolmAdjust <- function(p_values) {
    if (!is.numeric(p_values) || any(!is.finite(p_values)) ||
        any(p_values < 0 | p_values > 1)) {
        publicationStatisticsAbort("P-values are invalid")
    }
    stats::p.adjust(p_values, method = "holm")
}

publicationLogSd <- function(values) {
    values <- as.numeric(unlist(values, use.names = FALSE))
    if (length(values) < 2L || any(!is.finite(values)) || any(values <= 0)) {
        publicationStatisticsAbort("Precision source values are invalid")
    }
    stats::sd(log(values))
}

publicationPrecisionRequiredPairs <- function(
    log_sd,
    ratio_half_width = 1.10,
    confidence = 0.95
) {
    if (!publicationScalarNumber(log_sd) || log_sd < 0 ||
        !publicationScalarNumber(ratio_half_width, positive = TRUE) ||
        ratio_half_width <= 1) {
        publicationStatisticsAbort("Precision parameters are invalid")
    }
    z <- stats::qnorm(1 - (1 - confidence) / 2)
    ceiling((z * log_sd / log(ratio_half_width))^2)
}

publicationRoundPairCount <- function(required, minimum = 30L, maximum = 60L) {
    if (!publicationScalarNumber(required) || required < 0 ||
        minimum %% 6L != 0L || maximum %% 6L != 0L || minimum > maximum) {
        publicationStatisticsAbort("Pair-count bounds are invalid")
    }
    selected <- max(minimum, 6L * ceiling(required / 6L))
    list(
        required_pairs = as.integer(required),
        selected_pairs = if (selected <= maximum) as.integer(selected) else NULL,
        maximum_pairs = as.integer(maximum),
        status = if (selected <= maximum) "precision_satisfied" else
            "underpowered_at_maximum_non_promotional"
    )
}

publicationBaselinePrecisionEvidence <- function(paths) {
    records <- lapply(paths, \(path) {
        value <- publicationReadJson(path)
        summary <- value$summary
        list(
            path = path,
            sha256 = publicationFileDigest(path),
            workload_id = value$workload_id,
            repetitions = summary$completed,
            metrics = list(
                peak_charged_proxy = summary$peak_tree_rss_bytes,
                retained_charged_proxy = summary$retained_tree_rss_bytes,
                elapsed_seconds = summary$elapsed_seconds
            )
        )
    })
    records
}

publicationSelectPrecisionPairCount <- function(
    baseline_records,
    null_calibration,
    minimum = 30L,
    maximum = 60L,
    ratio_half_width = 1.10
) {
    metric_sds <- list()
    for (record in baseline_records) {
        for (metric in names(record$metrics)) {
            key <- paste(record$workload_id, metric, sep = "::")
            metric_sds[[key]] <- publicationLogSd(record$metrics[[metric]])
        }
    }
    null_metrics <- c(
        "peak_charged_memory_bytes",
        "retained_charged_memory_bytes",
        "elapsed_seconds",
        "cpu_usage_seconds"
    )
    for (metric in null_metrics) {
        values <- vapply(
            null_calibration$runs,
            \(run) as.numeric(run$metrics[[metric]]),
            numeric(1)
        )
        metric_sds[[paste0("null::", metric)]] <- publicationLogSd(values)
    }
    conservative_log_ratio_sd <- sqrt(2) * max(unlist(metric_sds), na.rm = TRUE)
    required <- publicationPrecisionRequiredPairs(
        conservative_log_ratio_sd,
        ratio_half_width
    )
    selection <- publicationRoundPairCount(required, minimum, maximum)
    list(
        schema = "multischolar.omics_publication_precision_successor",
        schema_version = "1.0.0",
        owner_ticket_id = "OMICS-ART-063",
        evidence_class = "immutable_pre_candidate_memory_and_null",
        candidate_loaded = FALSE,
        method = list(
            estimator = "conservative_independent_log_ratio_sd",
            formula = "sqrt(2) * max(sd(log(metric)))",
            confidence = 0.95,
            ratio_half_width = ratio_half_width,
            round_to_multiple = 6,
            minimum_pairs = minimum,
            maximum_pairs = maximum
        ),
        source_records = lapply(baseline_records, \(record) {
            record[c("path", "sha256", "workload_id", "repetitions")]
        }),
        null_calibration_id = null_calibration$calibration_id,
        null_calibration_digest = publicationObjectDigest(null_calibration),
        metric_log_sds = metric_sds,
        conservative_log_ratio_sd = conservative_log_ratio_sd,
        required_pairs = selection$required_pairs,
        selected_pairs = selection$selected_pairs,
        status = selection$status,
        promotion_authority = FALSE,
        successor_required_before_campaign = TRUE
    )
}

publicationNullGovernanceBindings <- function() {
    protocol_path <- "tests/testdata/omics-performance/protocol-v1.json"
    cache_path <- publicationCacheContractPath()
    host_path <- "tests/testdata/omics-performance/host-preflight-contract-v2.json"
    protocol <- publicationReadJson(protocol_path)
    cache <- publicationReadJson(cache_path)
    host_contract <- publicationReadJson(host_path)
    list(
        protocol = list(
            protocol_id = protocol$protocol_id,
            path = protocol_path,
            sha256 = publicationFileDigest(protocol_path)
        ),
        metric_contract = publicationMetricContractBinding(),
        cache_contract = list(
            cache_contract_id = cache$cache_contract_id,
            path = cache_path,
            sha256 = publicationFileDigest(cache_path)
        ),
        host_preflight_contract = list(
            preflight_contract_id = host_contract$preflight_contract_id,
            path = host_path,
            sha256 = publicationFileDigest(host_path)
        )
    )
}

publicationValidateNullCalibrationHeader <- function(record) {
    expected <- c(
        "schema", "schema_version", "calibration_id", "generated_at",
        "candidate_loaded", "implementation_sources", "governance_bindings",
        "method", "host", "runs", "status"
    )
    publicationRequireNames(record, expected, "Null cgroup calibration")
    publicationValidateKernelSourceBindings(record$implementation_sources)
    publicationRequireNames(
        record$governance_bindings,
        names(publicationNullGovernanceBindings()),
        "Calibration governance bindings"
    )
    method_fields <- c(
        "owned_cgroup_v2", "sampling_interval_ms",
        "disk_sampling_interval_ms", "retained_dwell_seconds",
        "retained_window_seconds", "retention_acknowledgement",
        "maximum_boundary_bracket_seconds",
        "retained_boundary_tolerance_ms",
        "candidate_loaded", "worker_thread_environment"
    )
    thread_names <- c(
        "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
        "ARROW_NUM_THREADS", "DUCKDB_THREADS"
    )
    publicationRequireNames(record$method, method_fields, "Calibration method")
    publicationRequireNames(
        record$method$worker_thread_environment,
        thread_names,
        "Calibration thread environment"
    )
    thread_values <- unlist(
        record$method$worker_thread_environment,
        use.names = FALSE
    )
    valid <- identical(
        record$schema,
        "multischolar.omics_publication_null_calibration"
    ) && identical(record$schema_version, "1.0.0") &&
        publicationScalarString(record$calibration_id) &&
        publicationScalarString(record$generated_at) &&
        identical(record$status, "passed") && !isTRUE(record$candidate_loaded) &&
        length(record$runs) == 30L &&
        identical(record$governance_bindings, publicationNullGovernanceBindings()) &&
        isTRUE(record$method$owned_cgroup_v2) &&
        !isTRUE(record$method$candidate_loaded) &&
        record$method$sampling_interval_ms == 20 &&
        record$method$disk_sampling_interval_ms == 500 &&
        record$method$retained_dwell_seconds == 5 &&
        record$method$retained_window_seconds == 2 &&
        identical(record$method$retention_acknowledgement, "fifo_v1") &&
        record$method$maximum_boundary_bracket_seconds == 0.5 &&
        record$method$retained_boundary_tolerance_ms == 100 &&
        all(is.finite(thread_values)) && all(thread_values == 1) &&
        identical(record$host$cgroup$version, 2L) &&
        isTRUE(record$host$cgroup$systemd_user_available) &&
        identical(
            record$host$cgroup$current_path,
            "<redacted-owned-user-cgroup>"
        )
    run_ids <- vapply(record$runs, `[[`, character(1), "run_id")
    valid <- valid && identical(run_ids, sprintf("null-%03d", seq_len(30L)))
    if (!valid) {
        publicationStatisticsAbort("Null cgroup calibration contract differs")
    }
    invisible(record)
}

publicationNullRunFields <- function() {
    c(
        "run_id", "status", "metric_contract", "exit_status", "timed_out",
        "safety_aborted", "safety_reason", "cgroup_observed", "cgroup_lost",
        "expected_root_pid", "pid_ownership_valid", "unexpected_pids",
        "indeterminate_pids", "retention_acknowledged", "unit_result",
        "worker_resource_ledger",
        "worker_resource_ledger_valid", "cleanup", "publication_certifiable",
        "retention_state", "retention_state_valid", "retained_window",
        "retained_diagnostics", "metrics_valid", "retained_sample_count",
        "metrics", "phase_evidence", "worker", "raw_measurement_sha256"
    )
}

publicationNullMetricFields <- function() {
    list(
        positive = c(
            "elapsed_seconds", "peak_charged_memory_bytes",
            "peak_sampled_charged_memory_bytes",
            "retained_charged_memory_bytes", "peak_anonymous_memory_bytes",
            "retained_anonymous_memory_bytes", "peak_pss_bytes",
            "retained_pss_bytes", "peak_uss_bytes", "retained_uss_bytes",
            "peak_rss_bytes", "retained_rss_bytes", "cpu_usage_seconds",
            "maximum_process_count", "peak_file_descriptors",
            "retained_file_descriptors", "peak_thread_tasks",
            "retained_thread_tasks", "peak_logical_disk_bytes",
            "peak_allocated_disk_bytes", "final_logical_disk_bytes",
            "final_allocated_disk_bytes", "final_file_count", "sample_count",
            "cgroup_lifecycle_elapsed_seconds", "phase_cpu_seconds",
            "phase_process_self_cpu_seconds",
            "primary_work_units_per_wall_second",
            "primary_work_units_per_cpu_second"
        ),
        nonnegative = c(
            "peak_file_memory_bytes", "retained_file_memory_bytes",
            "peak_kernel_memory_bytes", "retained_kernel_memory_bytes",
            "peak_slab_memory_bytes", "retained_slab_memory_bytes",
            "peak_sampled_accounted_memory_bytes",
            "retained_accounted_memory_bytes",
            "peak_reconciliation_absolute_residual_bytes",
            "retained_reconciliation_absolute_residual_bytes",
            "peak_swap_bytes", "cpu_user_seconds", "cpu_system_seconds",
            "peak_kernel_locks", "retained_kernel_locks",
            "harness_elapsed_seconds"
        ),
        retained = c(
            "sample_count", "observed_window_seconds",
            "charged_memory_slope_bytes_per_second",
            "anonymous_memory_slope_bytes_per_second", "cpu_activity_seconds",
            "io_read_activity_bytes", "io_write_activity_bytes",
            "background_activity_observed", "maximum_swap_bytes",
            "maximum_process_count", "trace_sha256"
        )
    )
}

publicationValidateNullRunTerminal <- function(run) {
    fields <- publicationNullMetricFields()
    publicationRequireNames(run, publicationNullRunFields(), "Null cgroup run")
    publicationRequireNames(
        run$retained_diagnostics,
        fields$retained,
        "Retained diagnostics"
    )
    publicationRequireNames(run$cleanup, c(
        "valid", "unit_active", "cgroup_exists", "elapsed_seconds"
    ), "Null cleanup")
    publicationRequireNames(run$retained_window, c(
        "valid", "expected_start_seconds", "expected_stop_seconds",
        "observed_start_seconds", "observed_stop_seconds",
        "pre_window_sample_seconds", "post_boundary_sample_seconds",
        "lower_boundary_bracketed", "lower_boundary_bracket_seconds",
        "upper_boundary_bracketed", "upper_boundary_bracket_seconds",
        "boundary_tolerance_seconds",
        "maximum_boundary_bracket_seconds"
    ), "Null retained window")
    raw_digest_valid <- publicationScalarString(run$raw_measurement_sha256) &&
        grepl("^[0-9a-f]{64}$", run$raw_measurement_sha256)
    trace_digest_valid <- publicationScalarString(
        run$retained_diagnostics$trace_sha256
    ) && grepl("^[0-9a-f]{64}$", run$retained_diagnostics$trace_sha256)
    terminal_valid <- identical(run$status, "passed") &&
        publicationScalarNumber(run$exit_status) && run$exit_status == 0 &&
        !isTRUE(run$timed_out) && !isTRUE(run$safety_aborted) &&
        is.null(run$safety_reason) && isTRUE(run$cgroup_observed) &&
        !isTRUE(run$cgroup_lost) &&
        publicationScalarNumber(run$expected_root_pid, positive = TRUE) &&
        isTRUE(run$pid_ownership_valid) && !length(run$unexpected_pids) &&
        isTRUE(run$retention_acknowledged) &&
        isTRUE(run$unit_result$available) &&
        identical(run$unit_result$result, "success") &&
        !isTRUE(run$unit_result$oom_killed) &&
        isTRUE(run$worker_resource_ledger_valid) &&
        publicationWorkerResourceLedgerValid(run$worker_resource_ledger) &&
        isTRUE(run$cleanup$valid) && !isTRUE(run$cleanup$unit_active) &&
        !isTRUE(run$cleanup$cgroup_exists) &&
        publicationScalarNumber(run$cleanup$elapsed_seconds) &&
        isTRUE(run$publication_certifiable) &&
        isTRUE(run$retention_state_valid) &&
        publicationRetentionStateValid(run$retention_state) &&
        isTRUE(run$retained_window$valid) &&
        isTRUE(run$retained_window$lower_boundary_bracketed) &&
        isTRUE(run$retained_window$upper_boundary_bracketed) &&
        publicationScalarNumber(
            run$retained_window$lower_boundary_bracket_seconds
        ) && publicationScalarNumber(
            run$retained_window$upper_boundary_bracket_seconds
        ) && run$retained_window$lower_boundary_bracket_seconds <=
            run$retained_window$maximum_boundary_bracket_seconds &&
        run$retained_window$upper_boundary_bracket_seconds <=
            run$retained_window$maximum_boundary_bracket_seconds &&
        isTRUE(run$metrics_valid) &&
        isTRUE(run$phase_evidence$valid) && raw_digest_valid &&
        trace_digest_valid &&
        identical(run$metric_contract, publicationMetricContractBinding())
    if (!terminal_valid) {
        publicationStatisticsAbort("Null cgroup run is not certifiable")
    }
    invisible(run)
}

publicationValidateNullRunMetrics <- function(run) {
    fields <- publicationNullMetricFields()
    positive_valid <- all(vapply(
        run$metrics[fields$positive],
        publicationScalarNumber,
        logical(1),
        positive = TRUE
    ))
    nonnegative_valid <- all(vapply(
        run$metrics[fields$nonnegative],
        \(value) publicationScalarNumber(value) && value >= 0,
        logical(1)
    ))
    io_valid <- all(vapply(c(run$metrics$io, run$metrics$phase_io), \(value) {
        publicationScalarNumber(value) && value >= 0
    }, logical(1)))
    event_valid <- all(vapply(run$metrics$memory_events, \(value) {
        publicationScalarNumber(value) && value >= 0
    }, logical(1))) && run$metrics$memory_events$oom == 0 &&
        run$metrics$memory_events$oom_kill == 0
    diagnostics <- run$retained_diagnostics
    diagnostics_valid <- all(vapply(diagnostics[c(
        "observed_window_seconds", "charged_memory_slope_bytes_per_second",
        "anonymous_memory_slope_bytes_per_second", "cpu_activity_seconds",
        "io_read_activity_bytes", "io_write_activity_bytes",
        "maximum_swap_bytes"
    )], publicationScalarNumber, logical(1))) &&
        diagnostics$sample_count == run$retained_sample_count &&
        !isTRUE(diagnostics$background_activity_observed) &&
        diagnostics$maximum_swap_bytes == 0 &&
        diagnostics$maximum_process_count == 1
    phase <- run$phase_evidence
    phase_valid <- publicationScalarNumber(
        phase$phase_cgroup_cpu_seconds,
        positive = TRUE
    ) && publicationScalarNumber(
        phase$phase_process_self_cpu_seconds,
        positive = TRUE
    ) && identical(phase$work_unit_id, run$worker$work_unit_id) &&
        phase$work_count == run$worker$work_count
    worker_valid <- identical(
        run$worker$schema,
        "multischolar.omics_publication_null_worker"
    ) && identical(run$worker$schema_version, "1.0.0") &&
        identical(run$worker$status, "passed") &&
        publicationScalarNumber(run$worker$checksum)
    if (!positive_valid || !nonnegative_valid || !io_valid ||
        !event_valid || !diagnostics_valid || !phase_valid || !worker_valid) {
        publicationStatisticsAbort("Null cgroup metrics are incomplete")
    }
    invisible(run)
}

publicationValidateNullCalibration <- function(record) {
    publicationValidateNullCalibrationHeader(record)
    for (run in record$runs) {
        publicationValidateNullRunTerminal(run)
        publicationValidateNullRunMetrics(run)
    }
    invisible(record)
}

publicationValidatePrecisionSuccessor <- function(record, null_calibration) {
    publicationValidateNullCalibration(null_calibration)
    required <- c(
        "schema", "schema_version", "owner_ticket_id", "evidence_class",
        "candidate_loaded", "method", "source_records", "null_calibration_id",
        "null_calibration_digest", "metric_log_sds",
        "conservative_log_ratio_sd", "required_pairs", "selected_pairs",
        "status", "promotion_authority", "successor_required_before_campaign",
        "precision_successor_id", "protocol_id", "null_calibration_path",
        "null_calibration_file_sha256"
    )
    publicationRequireNames(record, required, "Precision successor")
    publicationRequireNames(record$method, c(
        "estimator", "formula", "confidence", "ratio_half_width",
        "round_to_multiple", "minimum_pairs", "maximum_pairs"
    ), "Precision method")
    protocol <- publicationReadJson(
        "tests/testdata/omics-performance/protocol-v1.json"
    )
    common_invalid <- !identical(
        record$schema,
        "multischolar.omics_publication_precision_successor"
    ) || !identical(record$schema_version, "1.0.0") ||
        !identical(record$owner_ticket_id, "OMICS-ART-063") ||
        !identical(
            record$evidence_class,
            "immutable_pre_candidate_memory_and_null"
        ) || !publicationScalarString(record$precision_successor_id) ||
        !identical(record$protocol_id, protocol$protocol_id) ||
        !identical(record$null_calibration_id, null_calibration$calibration_id) ||
        !publicationScalarString(record$null_calibration_path) ||
        isTRUE(record$candidate_loaded) || isTRUE(record$promotion_authority) ||
        !isTRUE(record$successor_required_before_campaign) ||
        !identical(
            record$method$estimator,
            "conservative_independent_log_ratio_sd"
        ) || !identical(
            record$method$formula,
            "sqrt(2) * max(sd(log(metric)))"
        ) || record$method$confidence != 0.95 ||
        record$method$ratio_half_width != 1.10 ||
        record$method$minimum_pairs != 30L ||
        record$method$maximum_pairs != 60L ||
        record$method$round_to_multiple != 6L
    satisfied <- identical(record$status, "precision_satisfied") &&
        publicationScalarNumber(record$selected_pairs, positive = TRUE) &&
        record$selected_pairs >= record$method$minimum_pairs &&
        record$selected_pairs <= record$method$maximum_pairs &&
        record$selected_pairs %% record$method$round_to_multiple == 0L &&
        record$required_pairs <= record$method$maximum_pairs
    underpowered <- identical(
        record$status,
        "underpowered_at_maximum_non_promotional"
    ) && is.null(record$selected_pairs) &&
        record$required_pairs > record$method$maximum_pairs
    if (common_invalid || (!satisfied && !underpowered)) {
        publicationStatisticsAbort("Precision successor decision differs")
    }
    source_fields <- c("path", "sha256", "workload_id", "repetitions")
    source_paths <- vapply(record$source_records, \(source) {
        publicationRequireNames(source, source_fields, "Precision source")
        source$path
    }, character(1))
    metric_sds <- unlist(record$metric_log_sds, use.names = TRUE)
    if (length(source_paths) != 11L || anyDuplicated(source_paths) ||
        any(!is.finite(metric_sds)) || any(metric_sds < 0) ||
        !publicationScalarNumber(record$conservative_log_ratio_sd, positive = TRUE) ||
        !publicationScalarNumber(record$required_pairs, positive = TRUE)) {
        publicationStatisticsAbort("Precision source or variance differs")
    }
    if (!identical(
        publicationFileDigest(record$null_calibration_path),
        record$null_calibration_file_sha256
    ) || !identical(
        publicationObjectDigest(null_calibration),
        record$null_calibration_digest
    )) {
        publicationStatisticsAbort("Null calibration binding differs")
    }
    for (source in record$source_records) {
        if (!identical(publicationFileDigest(source$path), source$sha256) ||
            !publicationScalarString(source$workload_id) ||
            !publicationScalarNumber(source$repetitions, positive = TRUE)) {
            publicationStatisticsAbort("Precision source digest differs")
        }
    }
    baselines <- publicationBaselinePrecisionEvidence(vapply(
        record$source_records,
        `[[`,
        character(1),
        "path"
    ))
    replay <- publicationSelectPrecisionPairCount(
        baselines,
        null_calibration,
        minimum = record$method$minimum_pairs,
        maximum = record$method$maximum_pairs,
        ratio_half_width = record$method$ratio_half_width
    )
    if (!identical(replay$required_pairs, record$required_pairs) ||
        !identical(replay$selected_pairs, record$selected_pairs) ||
        !isTRUE(all.equal(
            replay$conservative_log_ratio_sd,
            record$conservative_log_ratio_sd,
            tolerance = 1e-12
        ))) {
        publicationStatisticsAbort("Precision successor does not replay")
    }
    invisible(record)
}
