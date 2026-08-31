publicationStatisticsFixture <- function(ratios = c(0.5, 0.75, 1.0)) {
    rows <- list()
    index <- 1L
    for (project_index in seq_along(ratios)) {
        for (pair_index in seq_len(10L)) {
            pair_id <- paste0("project-", project_index, "::pair-", pair_index)
            rows[[index]] <- data.frame(
                project_id = paste0("project-", project_index),
                pair_id = pair_id,
                arm = "A",
                status = "passed",
                elapsed_seconds = 100,
                stringsAsFactors = FALSE
            )
            index <- index + 1L
            rows[[index]] <- data.frame(
                project_id = paste0("project-", project_index),
                pair_id = pair_id,
                arm = "B",
                status = "passed",
                elapsed_seconds = 100 * ratios[[project_index]],
                stringsAsFactors = FALSE
            )
            index <- index + 1L
        }
    }
    do.call(rbind, rows)
}

test_that("paired and hierarchical effects match hand calculations", {
    records <- publicationStatisticsFixture()
    pairs <- publicationPairMetricTable(records, "elapsed_seconds")
    result <- publicationSummarizeMetric(
        records,
        "elapsed_seconds",
        resamples = 2000L,
        seed = 17L
    )

    expect_equal(nrow(pairs), 30L)
    expect_equal(
        result$summary$geometric_mean_ratio,
        exp(mean(log(c(0.5, 0.75, 1.0))))
    )
    expect_identical(result$hierarchical$resampling_unit, "project_then_process_pair")
    expect_identical(result$hierarchical$projects, 3L)
    expect_lte(result$hierarchical$lower, result$hierarchical$geometric_mean_ratio)
    expect_gte(result$hierarchical$upper, result$hierarchical$geometric_mean_ratio)
    expect_identical(result$analysis_seed, 17L)
    expect_equal(nrow(result$all_observations), nrow(records))
    expect_equal(nrow(result$paired_observations), 30L)
    expect_equal(result$failure_summary$failed, 0L)
    expect_match(result$source_result_digest, "^[0-9a-f]{64}$")
    expect_equal(
        result$project_results[["project-1"]]$bootstrap$lower,
        0.5
    )
    expect_equal(
        result$project_results[["project-1"]]$bootstrap$upper,
        0.5
    )
})

test_that("hierarchical bootstrap gives projects equal weight after failures", {
    records <- publicationStatisticsFixture()
    records <- records[!(records$project_id == "project-2" &
        !records$pair_id %in% paste0("project-2::pair-", 1:3)), ]
    records <- records[!(records$project_id == "project-3" &
        !records$pair_id %in% paste0("project-3::pair-", 1:3)), ]
    pairs <- publicationPairMetricTable(records, "elapsed_seconds")
    result <- publicationBootstrapHierarchical(
        pairs,
        resamples = 1000L,
        seed = 19L
    )

    expected <- exp(mean(log(c(0.5, 0.75, 1.0))))
    pooled <- exp(mean(pairs$log_ratio))
    expect_equal(result$geometric_mean_ratio, expected)
    expect_false(isTRUE(all.equal(result$geometric_mean_ratio, pooled)))
    expect_identical(
        result$project_weighting,
        "equal_project_weight_after_within_project_resampling"
    )
})

test_that("statistics reject missing arms nonfinite values and pseudoreplication", {
    records <- publicationStatisticsFixture()
    missing <- records[-1L, ]
    expect_error(
        publicationPairMetricTable(missing, "elapsed_seconds"),
        class = "multischolar_publication_statistics_error"
    )

    nonfinite <- records
    nonfinite$elapsed_seconds[[1L]] <- NA_real_
    expect_error(
        publicationPairMetricTable(nonfinite, "elapsed_seconds"),
        class = "multischolar_publication_statistics_error"
    )

    pairs <- publicationPairMetricTable(records, "elapsed_seconds")
    one_project <- pairs[pairs$project_id == "project-1", ]
    expect_error(
        publicationBootstrapHierarchical(one_project, resamples = 1000L),
        class = "multischolar_publication_statistics_error"
    )

    host_mismatch <- records
    host_mismatch$host_id <- "primary-host"
    host_mismatch$host_id[[2L]] <- "secondary-host"
    expect_error(
        publicationPairMetricTable(host_mismatch, "elapsed_seconds"),
        class = "multischolar_publication_statistics_error"
    )
})

test_that("failed pairs remain reported without imputation or replacement", {
    records <- publicationStatisticsFixture()
    failed_pair <- records$pair_id[[1L]]
    failed_row <- records$pair_id == failed_pair & records$arm == "B"
    records$status[failed_row] <- "worker_crash"
    records$elapsed_seconds[failed_row] <- NA_real_

    result <- publicationSummarizeMetric(
        records,
        "elapsed_seconds",
        resamples = 1000L,
        seed = 29L
    )

    expect_identical(result$excluded_pair_count, 1L)
    expect_identical(result$excluded_pair_ids, list(failed_pair))
    expect_identical(result$failure_summary$failed, 1L)
    expect_identical(result$failure_summary$failed_pairs, 1L)
    expect_equal(nrow(result$paired_observations), 29L)
    expect_false(failed_pair %in% result$paired_observations$pair_id)
})

test_that("Holm adjustment uses the declared secondary family", {
    expect_equal(
        publicationHolmAdjust(c(0.01, 0.02, 0.5)),
        stats::p.adjust(c(0.01, 0.02, 0.5), method = "holm")
    )
    expect_error(
        publicationHolmAdjust(c(0.1, NA_real_)),
        class = "multischolar_publication_statistics_error"
    )
})

publicationQueryFixture <- function(queries_per_run = 200L) {
    rows <- list()
    index <- 1L
    digest <- paste(rep("a", 64L), collapse = "")
    for (project in seq_len(3L)) {
        for (pair in seq_len(2L)) {
            for (arm in c("A", "B")) {
                run_id <- paste(project, pair, arm, sep = "-")
                ratio <- if (arm == "B") 0.8 else 1
                for (query in seq_len(queries_per_run)) {
                    rows[[index]] <- data.frame(
                        project_id = paste0("project-", project),
                        pair_id = paste0("project-", project, "::pair-", pair),
                        arm = arm,
                        run_id = run_id,
                        status = "passed",
                        query_id = paste0("query-", query),
                        latency_seconds = ratio * (1 + query / 10000),
                        query_suite_elapsed_seconds = queries_per_run * 1.1,
                        returned_rows = 10,
                        returned_bytes = 1000,
                        output_digest = digest,
                        stringsAsFactors = FALSE
                    )
                    index <- index + 1L
                }
            }
        }
    }
    do.call(rbind, rows)
}

test_that("query statistics resample complete process runs", {
    queries <- publicationQueryFixture()
    result <- publicationSummarizeClusteredQueries(
        queries,
        resamples = 1000L,
        seed = 23L
    )

    expect_identical(result$resampling_unit, "fresh_process_run_pair")
    expect_false(result$within_run_queries_are_independent_replicates)
    expect_equal(nrow(result$run_observations), 12L)
    expect_equal(
        result$metric_results$query_p95_seconds$hierarchical$pairs,
        6L
    )
    expect_equal(
        result$metric_results$query_p95_seconds$hierarchical$projects,
        3L
    )
    expect_lt(
        result$metric_results$query_p95_seconds$hierarchical$geometric_mean_ratio,
        1
    )

    malformed <- queries
    malformed$output_digest[[1L]] <- "not-a-sha256"
    expect_error(
        publicationQueryRunTable(malformed),
        class = "multischolar_publication_statistics_error"
    )

    short <- publicationQueryFixture(199L)
    expect_true(all(publicationQueryRunTable(short)$status ==
        "query_contract_failed"))
})

test_that("precision selector is pre-candidate bounded and deterministic", {
    values <- c(100, 101, 99, 100.5, 99.5)
    baseline <- list(list(
        path = "baseline.json",
        sha256 = paste(rep("a", 64L), collapse = ""),
        workload_id = "baseline-workload",
        repetitions = 5L,
        metrics = list(
            peak_charged_proxy = values,
            retained_charged_proxy = values,
            elapsed_seconds = values
        )
    ))
    null <- list(
        calibration_id = "null-v1",
        runs = lapply(seq_len(30L), \(index) {
            list(metrics = list(
                peak_charged_memory_bytes = 100 + (index %% 3),
                retained_charged_memory_bytes = 99 + (index %% 3),
                elapsed_seconds = 1 + (index %% 3) / 100,
                cpu_usage_seconds = 0.5 + (index %% 3) / 100
            ))
        })
    )
    selected <- publicationSelectPrecisionPairCount(baseline, null)

    expect_identical(selected$status, "precision_satisfied")
    expect_identical(selected$selected_pairs, 30L)
    expect_false(selected$candidate_loaded)

    high_variance <- publicationRoundPairCount(61L)
    expect_identical(high_variance$status, "underpowered_at_maximum_non_promotional")
    expect_null(high_variance$selected_pairs)
})

test_that("real null calibration and precision successor replay exactly", {
    null <- publicationReadJson(
        "tests/testdata/omics-performance/null-cgroup-calibration-v6.json"
    )
    successor <- publicationReadJson(
        "tests/testdata/omics-performance/precision-successor-v6.json"
    )

    expect_silent(publicationValidateNullCalibration(null))
    expect_silent(publicationValidatePrecisionSuccessor(successor, null))
    expect_identical(successor$required_pairs, 31L)
    expect_identical(successor$selected_pairs, 36L)
    expect_length(successor$source_records, 11L)
    expect_false(successor$candidate_loaded)
    expect_false(successor$promotion_authority)

    drift <- publicationGovernanceCopy(successor)
    drift$selected_pairs <- 42L
    expect_error(
        publicationValidatePrecisionSuccessor(drift, null),
        class = "multischolar_publication_statistics_error"
    )

    candidate_leak <- publicationGovernanceCopy(successor)
    candidate_leak$candidate_loaded <- TRUE
    expect_error(
        publicationValidatePrecisionSuccessor(candidate_leak, null),
        class = "multischolar_publication_statistics_error"
    )
})

test_that("governance manifest v3 binds precision records and kernel sources", {
    manifest <- publicationReadJson(
        "tests/testdata/omics-performance/governance-manifest-v3.json"
    )

    expect_silent(publicationValidateGovernanceManifestV3(manifest))
    expect_identical(manifest$decisions$exact_primary_pair_count, 36L)
    expect_false(manifest$decisions$publication_host_certified)
    expect_false(manifest$decisions$campaign_execution_authorized)
    expect_false(manifest$decisions$promotion_authorized)

    drift <- publicationGovernanceCopy(manifest)
    drift$records[[1L]]$sha256 <- paste(rep("0", 64L), collapse = "")
    expect_error(
        publicationValidateGovernanceManifestV3(drift),
        class = "multischolar_publication_binding_error"
    )
})
