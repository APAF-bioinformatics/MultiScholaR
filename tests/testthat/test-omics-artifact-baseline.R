.omicsParityNamed <- function(values, field) {
    stats::setNames(values, vapply(values, `[[`, character(1), field))
}

test_that("DIA-NN parity manifest and memory oracle are versioned and exact", {
    manifest <- omicsParityReadManifest()
    oracle <- omicsParityReadOracle()

    expect_identical(manifest$schema_version, .OMICS_PARITY_MANIFEST_SCHEMA)
    expect_identical(oracle$schema_version, .OMICS_PARITY_ORACLE_SCHEMA)
    expect_identical(manifest$backend, "memory")
    expect_identical(oracle$backend, "memory")
    expect_identical(manifest$capability_id, "proteomics.diann.peptide.dia.v1")
    expect_match(oracle$oracle_origin, "pre-artifact", fixed = TRUE)

    scenario_ids <- vapply(manifest$scenarios, `[[`, character(1), "scenario_id")
    oracle_ids <- vapply(oracle$scenarios, `[[`, character(1), "scenario_id")
    expect_setequal(
        oracle_ids,
        c("dia_canonical_import", "dia_limpa_mnar_import")
    )
    expect_true(all(oracle_ids %in% scenario_ids))
    expect_identical(anyDuplicated(scenario_ids), 0L)
})

test_that("manifest separates canonical, MNAR, scaling, and private evidence", {
    manifest <- omicsParityReadManifest()
    scenarios <- .omicsParityNamed(manifest$scenarios, "scenario_id")

    canonical <- scenarios[["dia_canonical_import"]]
    limpa <- scenarios[["dia_limpa_mnar_import"]]
    scaling <- scenarios[["dia_generated_scaling_250x"]]
    private <- scenarios[["dia_private_representative"]]

    expect_identical(canonical$kind, "committed_fixture")
    expect_identical(limpa$kind, "committed_fixture")
    expect_identical(scaling$kind, "generated_scaling")
    expect_identical(private$kind, "optional_private")
    expect_identical(scaling$row_multiplier, 250L)
    expect_identical(private$fixture_env, "MULTISCHOLAR_DIA_BASELINE_FIXTURE")
    expect_null(private$fixture_path)
    expect_false(isTRUE(private$required))
    expect_match(private$privacy, "never persist source path", fixed = TRUE)

    for (scenario in list(canonical, limpa)) {
        expect_true(file.exists(omicsParityRepoPath(scenario$fixture_path)))
        expect_true(file.exists(omicsParityRepoPath(scenario$design_path)))
        expect_identical(scenario$parameters$backend, "memory")
    }

    canonical_data <- utils::read.delim(
        omicsParityRepoPath(canonical$fixture_path),
        check.names = FALSE
    )
    limpa_data <- utils::read.delim(
        omicsParityRepoPath(limpa$fixture_path),
        check.names = FALSE
    )
    expect_equal(nrow(canonical_data), 60L)
    expect_equal(sum(is.na(canonical_data$Precursor.Normalised)), 0L)
    expect_equal(nrow(limpa_data), 72L)
    expect_equal(sum(is.na(limpa_data$Precursor.Normalised)), 20L)
    expect_setequal(unique(limpa_data$Protein.Group), sprintf("P%05d", 1:6))
})

test_that("comparison contract names keys, ordering, tolerances, and special values", {
    contract <- omicsParityReadManifest()$comparison_contract

    expect_identical(
        unlist(contract$stable_keys, use.names = FALSE),
        c("Run", "Protein.Group", "Stripped.Sequence", "Precursor.Id")
    )
    expect_length(contract$ordered_columns, 16L)
    expect_true(all(c(
        "input_sha256", "output_sha256", "column_names", "column_classes",
        "na_quantity", "nan_quantity", "positive_infinite_quantity",
        "negative_infinite_quantity"
    ) %in% unlist(contract$exact_fields, use.names = FALSE)))
    expect_setequal(
        names(contract$numeric_fields),
        c("quantity_sum", "quantity_min", "quantity_max")
    )
    expect_true(all(vapply(contract$numeric_fields, \(tolerance) {
        tolerance$absolute >= 0 && tolerance$relative >= 0
    }, logical(1))))
    expect_setequal(
        names(contract$special_values),
        c("missing", "not_a_number", "positive_infinity", "negative_infinity", "signed_zero")
    )
})

test_that("pre-artifact DIA-NN imports reproduce the committed memory oracle", {
    manifest <- omicsParityReadManifest()
    oracle <- omicsParityReadOracle()
    scenarios <- .omicsParityNamed(manifest$scenarios, "scenario_id")
    expected <- .omicsParityNamed(oracle$scenarios, "scenario_id")

    for (scenario_id in names(expected)) {
        scenario <- scenarios[[scenario_id]]
        fixture_path <- omicsParityRepoPath(scenario$fixture_path)
        imported <- suppressMessages(importDIANNData(
            fixture_path,
            use_precursor_norm = isTRUE(scenario$parameters$use_precursor_norm)
        ))
        actual <- omicsParitySummarizeDiann(
            imported,
            fixture_path,
            manifest$comparison_contract
        )
        comparison <- omicsParityCompareSummary(
            actual,
            expected[[scenario_id]]$expected,
            manifest$comparison_contract
        )
        expect_true(
            comparison$equal,
            info = paste(unlist(comparison$errors, use.names = FALSE), collapse = "; ")
        )
    }
})

test_that("canonical table hashes are order-independent but value-sensitive", {
    manifest <- omicsParityReadManifest()
    contract <- manifest$comparison_contract
    fixture_path <- omicsParityRepoPath(
        manifest$scenarios[[1L]]$fixture_path
    )
    imported <- suppressMessages(importDIANNData(fixture_path))
    data <- as.data.frame(imported$data)
    stable_keys <- unlist(contract$stable_keys, use.names = FALSE)
    columns <- unlist(contract$ordered_columns, use.names = FALSE)
    expected_hash <- omicsParityTableSha256(data, stable_keys, columns)

    set.seed(104729L)
    shuffled <- data[sample.int(nrow(data)), , drop = FALSE]
    expect_identical(
        omicsParityTableSha256(shuffled, stable_keys, columns),
        expected_hash
    )
    changed <- data
    changed$Precursor.Normalised[[1L]] <- changed$Precursor.Normalised[[1L]] + 1
    expect_false(identical(
        omicsParityTableSha256(changed, stable_keys, columns),
        expected_hash
    ))
})

test_that("scientific corpus names separate methods and all downstream evidence", {
    manifest <- omicsParityReadManifest()
    targets <- .omicsParityNamed(manifest$scientific_targets, "target_id")

    expect_setequal(names(targets), c(
        "dia_browser_canonical", "dia_browser_limpa",
        "technical_replicate_peptide_imputation", "iq_and_limpa_rollup",
        "normalization_and_ruv", "differential_abundance_golden_master",
        "enrichment_doubles", "session_summary_and_reports"
    ))
    expect_false(identical(
        targets$technical_replicate_peptide_imputation$branch,
        targets$iq_and_limpa_rollup$branch
    ))
    expect_setequal(
        unlist(targets$iq_and_limpa_rollup$parameters$methods, use.names = FALSE),
        c("iq_maxlfq", "limpa_dpc_quant")
    )
    expect_setequal(
        vapply(manifest$cross_omic_regression_targets, `[[`, character(1), "target_id"),
        c("lfq_tmt", "metabolomics", "lipidomics")
    )
    for (target in manifest$scientific_targets) {
        expect_identical(target$parameters$backend, "memory")
        expect_true(nzchar(target$test_filter))
    }
})

test_that("performance gates predate dual-write and distinguish release metrics", {
    manifest <- omicsParityReadManifest()
    gates <- manifest$release_gates

    expect_gte(gates$minimum_peak_rss_reduction_fraction, 0.30)
    expect_gte(gates$minimum_retained_rss_reduction_fraction, 0.40)
    expect_lte(gates$maximum_runtime_ratio, 1.25)
    expect_lte(gates$maximum_bounded_query_p95_seconds, 0.50)
    expect_gte(gates$minimum_representative_repetitions, 5L)
    expect_match(gates$gate_policy, "representative-project", fixed = TRUE)
    expect_setequal(
        vapply(manifest$disk_categories, `[[`, character(1), "category"),
        c("diagnostics", "committed", "staging_snapshot", "duckdb_spill", "final", "harness")
    )
    expect_identical(manifest$execution$repetitions, 2L)
    expect_identical(
        unlist(manifest$execution$cache_sequence, use.names = FALSE),
        c("cold_uncontrolled", "warm_candidate")
    )

    configured_threads <- unlist(
        manifest$thread_controls$configured_environment,
        use.names = TRUE
    )
    expect_setequal(names(configured_threads), c(
        "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
        "ARROW_NUM_THREADS", "DUCKDB_THREADS"
    ))
    expect_true(all(configured_threads == 1L))
    iq_control <- manifest$thread_controls$method_level_limitations[[1L]]
    expect_identical(iq_control$component, "iq::process_long_format")
    expect_identical(iq_control$status, "reported_not_controlled")
    expect_match(iq_control$reason, "stdout", fixed = TRUE)
    expect_match(manifest$thread_controls$claim_policy, "reported separately", fixed = TRUE)
})

test_that("parity and result schemas reject unsupported future versions", {
    future_manifest <- jsonlite::read_json(omicsParityManifestPath(), simplifyVector = FALSE)
    future_manifest$schema_version <- "2.0.0"
    manifest_path <- tempfile(fileext = ".json")
    jsonlite::write_json(future_manifest, manifest_path, auto_unbox = TRUE)
    expect_error(
        omicsParityReadManifest(manifest_path),
        "Unsupported omics parity manifest schema version '2.0.0'",
        fixed = TRUE
    )

    future_result <- list(schema_version = "2.0.0")
    result_path <- tempfile(fileext = ".json")
    jsonlite::write_json(future_result, result_path, auto_unbox = TRUE)
    expect_error(
        omicsParityValidateResult(result_path),
        "Unsupported baseline result schema version '2.0.0'",
        fixed = TRUE
    )
})

test_that("fresh-process harness records RSS, disk, allocation, query, and native metrics", {
    skip_on_cran()
    skip_if_not_installed("processx")
    skip_if_not_installed("ps")

    output_dir <- tempfile("omics-artifact-baseline-")
    runner <- omicsParityRepoPath("tools", "profiling", "run_omics_artifact_baseline.R")
    process <- processx::run(
        file.path(R.home("bin"), "Rscript"),
        c(
            "--vanilla", runner,
            "--scenario", "dia_canonical_import",
            "--repetitions", "1",
            "--output-dir", output_dir
        ),
        wd = tempdir(),
        timeout = 180000,
        error_on_status = FALSE,
        echo = FALSE
    )
    expect_identical(process$status, 0L, info = process$stderr)

    result <- omicsParityValidateResult(file.path(output_dir, "baseline-result.json"))
    expect_length(result$runs, 1L)
    run <- result$runs[[1L]]
    expect_identical(run$status, "passed")
    expect_identical(run$worker$status, "passed")
    expect_true(run$worker$oracle_comparison$equal)
    expect_gt(run$metrics$peak_tree_rss_bytes, 0)
    expect_gt(run$metrics$retained_tree_rss_bytes, 0)
    expect_gt(run$metrics$sample_count, 1L)
    expect_false(run$worker$allocation_diagnostics$available)
    expect_match(
        run$worker$allocation_diagnostics$reason,
        "separate diagnostic process",
        fixed = TRUE
    )
    expect_length(run$worker$stages[[2L]]$results, 1L)
    expect_gt(run$worker$stages[[2L]]$results[[1L]]$rows, 0L)
    expect_gt(run$worker$stages[[2L]]$results[[1L]]$p95_seconds, 0)
    expect_lte(
        run$worker$stages[[2L]]$results[[1L]]$p95_seconds,
        result$release_gates$maximum_bounded_query_p95_seconds
    )
    expect_true(all(unlist(run$worker$thread_environment, use.names = FALSE) == "1"))
    expect_true(is.logical(run$worker$native_resources$arrow$available))
    expect_length(result$diagnostics, 1L)
    diagnostic <- result$diagnostics[[1L]]
    expect_identical(diagnostic$measurement_class, "allocation_copy_diagnostic")
    expect_true(diagnostic$worker$allocation_diagnostics$available)
    expect_gt(diagnostic$worker$allocation_diagnostics$allocation_records, 0L)
    expect_gt(diagnostic$metrics$peak_disk_category_bytes$diagnostics, 0)
    expect_identical(result$summary$completed, 1L)
    expect_identical(result$determinism[[1L]]$deterministic, TRUE)
    expect_identical(result$scenario_summaries$dia_canonical_import$completed, 1L)
    expect_identical(
        result$thread_controls$method_level_limitations[[1L]]$status,
        "reported_not_controlled"
    )
    expect_lt(
        run$metrics$peak_artifact_disk_bytes,
        diagnostic$metrics$peak_disk_bytes
    )
    expect_false(result$sampling$object_size_is_release_metric)
    expect_false(result$sampling$garbage_collection_is_release_proof)
    expect_true(result$sampling$rprofmem_is_allocation_diagnostic_only)
    expect_true(result$sampling$tracemem_is_copy_diagnostic_only)
})

test_that("private fixture evidence excludes source paths and biological identifiers", {
    skip_on_cran()
    skip_if_not_installed("processx")
    skip_if_not_installed("ps")

    manifest <- omicsParityReadManifest()
    canonical <- .omicsParityNamed(manifest$scenarios, "scenario_id")[[
        "dia_canonical_import"
    ]]
    private_dir <- tempfile("private-dia-source-marker-")
    dir.create(private_dir, recursive = TRUE)
    private_path <- file.path(private_dir, "private-dia-source-marker.tsv")
    expect_true(file.copy(omicsParityRepoPath(canonical$fixture_path), private_path))

    output_dir <- tempfile("omics-artifact-private-baseline-")
    runner <- omicsParityRepoPath("tools", "profiling", "run_omics_artifact_baseline.R")
    process <- processx::run(
        file.path(R.home("bin"), "Rscript"),
        c(
            "--vanilla", runner,
            "--scenario", "dia_private_representative",
            "--include-private", "true",
            "--diagnostics", "false",
            "--repetitions", "1",
            "--output-dir", output_dir
        ),
        wd = tempdir(),
        env = c(
            MULTISCHOLAR_DIA_BASELINE_FIXTURE = private_path,
            MULTISCHOLAR_BASELINE_FINGERPRINT_SALT = "test-only-private-salt"
        ),
        timeout = 180000,
        error_on_status = FALSE,
        echo = FALSE
    )
    expect_identical(process$status, 0L, info = process$stderr)

    result <- omicsParityValidateResult(file.path(output_dir, "baseline-result.json"))
    run <- result$runs[[1L]]
    summary <- run$worker$observed_summary
    expect_false(run$output_retained)
    expect_false(run$worker$fixture$source_path_retained)
    expect_identical(summary$run_count, 6L)
    expect_identical(summary$protein_count, 5L)
    expect_identical(summary$peptide_count, 10L)
    expect_null(summary$runs)
    expect_null(summary$proteins)
    expect_null(summary$peptides)
    expect_false(identical(summary$input_sha256, omicsParitySha256File(private_path)))

    persisted_paths <- list.files(output_dir, recursive = TRUE, full.names = TRUE)
    persisted_text <- unlist(lapply(persisted_paths, \(path) {
        paste(readLines(path, warn = FALSE), collapse = "\n")
    }), use.names = FALSE)
    expect_false(any(grepl(private_path, persisted_text, fixed = TRUE)))
    expect_false(any(grepl("WT_1|P00001|PEP1AK", persisted_text)))
})
