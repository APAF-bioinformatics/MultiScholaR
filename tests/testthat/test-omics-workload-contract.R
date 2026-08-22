.omicsWorkloadRepoPath <- function(...) {
    file.path(normalizePath(testthat::test_path("..", ".."), mustWork = TRUE), ...)
}

source(
    .omicsWorkloadRepoPath("tools", "profiling", "omics_workload_contract.R"),
    local = TRUE
)
source(
    .omicsWorkloadRepoPath("tools", "profiling", "run_omics_memory_baseline.R"),
    local = TRUE
)

.omicsWorkloadWriteAdapter <- function(path) {
    writeLines(c(
        "omicsWorkloadAdapter <- function() {",
        "    prepare <- function(context) {",
        "        contract <- context$contract",
        "        staging <- file.path(context$run_dir, 'staging')",
        "        dir.create(staging, recursive = TRUE, showWarnings = FALSE)",
        "        value_count <- contract$dimensions$feature_count *",
        "            contract$dimensions$sample_count",
        "        values <- stats::runif(value_count)",
        "        payload <- matrix(values, nrow = contract$dimensions$feature_count)",
        "        payload_path <- file.path(staging, 'payload.tsv')",
        "        truth_path <- file.path(staging, 'truth.json')",
        "        utils::write.table(",
        "            payload, payload_path, sep = '\\t', row.names = FALSE,",
        "            col.names = FALSE, quote = FALSE",
        "        )",
        "        truth <- list(",
        "            mean = mean(payload), rows = nrow(payload),",
        "            columns = ncol(payload)",
        "        )",
        "        jsonlite::write_json(",
        "            truth, truth_path, auto_unbox = TRUE, digits = 17",
        "        )",
        "        list(",
        "            payload_path = payload_path, truth_path = truth_path,",
        "            metadata = list(generated = TRUE)",
        "        )",
        "    }",
        "    run <- function(context) {",
        "        payload <- utils::read.delim(",
        "            context$payload_path, header = FALSE, check.names = FALSE",
        "        )",
        "        payload <- as.matrix(payload)",
        "        truth <- jsonlite::read_json(context$truth_path, simplifyVector = TRUE)",
        "        started <- proc.time()[['elapsed']]",
        "        selected <- payload[seq_len(min(2L, nrow(payload))), , drop = FALSE]",
        "        query_seconds <- proc.time()[['elapsed']] - started + .Machine$double.eps",
        "        list(",
        "            status = 'passed',",
        "            workflow_evidence = list(",
        "                truth_valid = isTRUE(all.equal(",
        "                    mean(payload), truth$mean, tolerance = 1e-15",
        "                )),",
        "                rows = nrow(payload), columns = ncol(payload)",
        "            ),",
        "            query_evidence = list(list(",
        "                query_id = 'head', rows = nrow(selected),",
        "                p95_seconds = query_seconds",
        "            )),",
        "            session_evidence = list(status = 'closed', resource_count = 0L),",
        "            report_evidence = list(status = 'not_required', file_count = 0L),",
        "            retained = payload",
        "        )",
        "    }",
        "    list(",
        "        adapter_id = 'test.synthetic.v1',",
        "        adapter_version = '1.0.0',",
        "        supported_omics = 'test', prepare = prepare, run = run",
        "    )",
        "}"
    ), path)
    path
}

.omicsWorkloadContract <- function(adapter_path, seed = 104729L) {
    list(
        schema = .OMICS_WORKLOAD_SCHEMA,
        schema_version = .OMICS_WORKLOAD_SCHEMA_VERSION,
        workload_id = "test.synthetic.small.v1",
        capability = list(
            capability_id = "test.synthetic.matrix.v1",
            omic_type = "test",
            input_format = "synthetic",
            data_level = "feature",
            acquisition_mode = "simulated"
        ),
        adapter = list(
            adapter_id = "test.synthetic.v1",
            adapter_version = "1.0.0",
            source_sha256 = omicsWorkloadFileDigest(adapter_path)
        ),
        rng = list(
            kind = as.list(c("Mersenne-Twister", "Inversion", "Rejection")),
            seed = seed
        ),
        dimensions = list(feature_count = 4L, sample_count = 3L, assay_count = 1L),
        assay_mix = list(test = 1L),
        distributions = list(intensity = "uniform_test_only"),
        missingness = list(mechanism = "none"),
        censoring = list(mechanism = "none"),
        duplicate_policy = list(mode = "none"),
        sample_perturbations = list(mode = "none"),
        adapter_parameters = list(),
        privacy = list(
            classification = "public_synthetic",
            source = "deterministic_test_adapter",
            scale_metadata = NULL
        ),
        execution = list(
            repetitions = 5L,
            sampling_interval_ms = 10L,
            retention_window_ms = 250L,
            timeout_ms = 5000L,
            maximum_retained_samples = 5L,
            threads = 1L,
            locale = "C",
            timezone = "UTC",
            cache_state = "cold_uncontrolled"
        ),
        expected_digests = list(
            payload_sha256 = paste(rep("0", 64L), collapse = ""),
            truth_sha256 = paste(rep("0", 64L), collapse = "")
        )
    )
}

.omicsWorkloadFreezeDigests <- function(contract, adapter_path) {
    adapter_environment <- new.env(parent = baseenv())
    sys.source(adapter_path, envir = adapter_environment)
    adapter <- adapter_environment$omicsWorkloadAdapter()
    run_dir <- tempfile("omics-workload-bootstrap-")
    dir.create(run_dir, recursive = TRUE)
    do.call(RNGkind, as.list(unlist(contract$rng$kind, use.names = FALSE)))
    set.seed(as.integer(contract$rng$seed))
    prepared <- adapter$prepare(list(
        contract = contract,
        run_dir = run_dir,
        repetition = 1L
    ))
    contract$expected_digests <- list(
        payload_sha256 = omicsWorkloadFileDigest(prepared$payload_path),
        truth_sha256 = omicsWorkloadFileDigest(prepared$truth_path)
    )
    contract
}

test_that("workload contracts reject schema, adapter, and privacy drift", {
    adapter_path <- .omicsWorkloadWriteAdapter(tempfile(fileext = ".R"))
    contract <- .omicsWorkloadFreezeDigests(
        .omicsWorkloadContract(adapter_path),
        adapter_path
    )
    expect_invisible(omicsWorkloadValidateContract(contract))

    future <- contract
    future$schema_version <- "2.0.0"
    expect_error(
        omicsWorkloadValidateContract(future),
        class = "omics_workload_contract_error"
    )
    unknown <- contract
    unknown$undeclared <- TRUE
    expect_error(
        omicsWorkloadValidateContract(unknown),
        class = "omics_workload_contract_error"
    )
    leaked <- contract
    leaked$adapter_parameters$source_path <- "/private/input.tsv"
    expect_error(
        omicsWorkloadValidateContract(leaked),
        class = "omics_workload_privacy_error"
    )
    empirical <- contract
    empirical$distributions$empirical_quantiles <- c(1, 2, 3)
    expect_error(
        omicsWorkloadValidateContract(empirical),
        class = "omics_workload_privacy_error"
    )
    changed_adapter <- contract
    changed_adapter$adapter$adapter_version <- "2.0.0"
    expect_error(
        omicsWorkloadLoadAdapter(adapter_path, changed_adapter),
        class = "omics_workload_binding_error"
    )

    payload_path <- tempfile(fileext = ".tsv")
    truth_path <- tempfile(fileext = ".json")
    writeLines("payload", payload_path)
    writeLines("truth", truth_path)
    mismatched <- contract
    mismatched$expected_digests$payload_sha256 <- paste(
        rep("f", 64L),
        collapse = ""
    )
    expect_error(
        omicsWorkloadEvidenceBinding(
            contract = mismatched,
            adapter_path = adapter_path,
            payload_path = payload_path,
            truth_path = truth_path,
            code_revision = list(git_sha = "test"),
            environment = list(r_version = R.version.string)
        ),
        class = "omics_workload_binding_error"
    )
})

test_that("contract and generated evidence separate changed seeds", {
    adapter_path <- .omicsWorkloadWriteAdapter(tempfile(fileext = ".R"))
    first <- .omicsWorkloadFreezeDigests(
        .omicsWorkloadContract(adapter_path, seed = 104729L),
        adapter_path
    )
    repeated <- .omicsWorkloadFreezeDigests(
        .omicsWorkloadContract(adapter_path, seed = 104729L),
        adapter_path
    )
    changed <- .omicsWorkloadFreezeDigests(
        .omicsWorkloadContract(adapter_path, seed = 104731L),
        adapter_path
    )

    expect_identical(first$expected_digests, repeated$expected_digests)
    expect_identical(omicsWorkloadDigest(first), omicsWorkloadDigest(repeated))
    expect_false(identical(first$expected_digests, changed$expected_digests))
    expect_false(identical(omicsWorkloadDigest(first), omicsWorkloadDigest(changed)))
})

test_that("shared memory runner measures five deterministic fresh processes", {
    skip_on_cran()
    skip_if_not_installed("processx")
    skip_if_not_installed("ps")

    adapter_path <- .omicsWorkloadWriteAdapter(tempfile(fileext = ".R"))
    contract <- .omicsWorkloadFreezeDigests(
        .omicsWorkloadContract(adapter_path),
        adapter_path
    )
    contract_path <- tempfile(fileext = ".json")
    jsonlite::write_json(
        contract,
        contract_path,
        auto_unbox = TRUE,
        pretty = TRUE,
        null = "null",
        digits = 17
    )
    output_dir <- tempfile("omics-memory-baseline-")
    process <- processx::run(
        file.path(R.home("bin"), "Rscript"),
        c(
            "--vanilla",
            .omicsWorkloadRepoPath(
                "tools",
                "profiling",
                "run_omics_memory_baseline.R"
            ),
            "--contract",
            contract_path,
            "--adapter",
            adapter_path,
            "--output-dir",
            output_dir,
            "--install-package",
            "false",
            "--diagnostics",
            "true"
        ),
        wd = .omicsWorkloadRepoPath(),
        timeout = 120000,
        error_on_status = FALSE,
        echo = FALSE
    )
    expect_identical(process$status, 0L, info = process$stderr)

    result <- jsonlite::read_json(
        file.path(output_dir, "baseline-result.json"),
        simplifyVector = FALSE
    )
    failed_runs <- Filter(\(run) !identical(run$status, "passed"), result$runs)
    expect_identical(
        result$status,
        "passed",
        info = jsonlite::toJSON(failed_runs, auto_unbox = TRUE)
    )
    expect_length(result$runs, 5L)
    expect_length(result$diagnostics, 1L)
    expect_true(all(unlist(result$determinism, use.names = FALSE)))
    expect_identical(result$preparation$fresh_processes, 2L)
    expect_true(result$preparation$deterministic)
    run_checks <- vapply(result$runs, \(run) {
        identical(run$status, "passed") &&
            run$metrics$peak_tree_rss_bytes > 0 &&
            run$metrics$retained_tree_rss_bytes > 0 &&
            run$metrics$committed_input_bytes > 0 &&
            run$metrics$committed_input_file_count == 2L &&
            all(c(
                "committed",
                "staging_snapshot",
                "duckdb_spill"
            ) %in% names(run$metrics$peak_disk_category_bytes)) &&
            length(run$samples) <= 5L &&
            isTRUE(run$preparation$outside_measured_interval) &&
            isTRUE(run$preparation$fresh_process) &&
            !isTRUE(run$worker$native_resources$arrow$loaded) &&
            !isTRUE(run$worker$native_resources$duckdb$loaded) &&
            run$worker$query_evidence[[1L]]$p95_seconds >= 0 &&
            run$worker$session_evidence$resource_count == 0L &&
            run$worker$report_evidence$file_count == 0L &&
            !length(run$validation$undeclared_files)
    }, logical(1))
    expect_true(
        all(run_checks),
        info = jsonlite::toJSON(result$runs[!run_checks], auto_unbox = TRUE)
    )
    expect_true(result$diagnostics[[1L]]$worker$allocation_diagnostics$available)
})

test_that("worker and file evidence fail closed when incomplete", {
    binding <- list(
        binding_sha256 = paste(rep("a", 64L), collapse = ""),
        payload_sha256 = paste(rep("b", 64L), collapse = ""),
        truth_sha256 = paste(rep("c", 64L), collapse = "")
    )
    incomplete <- list(status = "passed")
    expect_false(omicsMemoryValidateWorker(incomplete, binding)$valid)

    mixed_binding <- list(
        schema = .OMICS_MEMORY_WORKER_SCHEMA,
        schema_version = .OMICS_MEMORY_WORKER_VERSION,
        status = "passed",
        binding_sha256 = paste(rep("d", 64L), collapse = ""),
        payload_sha256 = binding$payload_sha256,
        truth_sha256 = binding$truth_sha256,
        retained_rss_bytes = 1
    )
    expect_false(omicsMemoryValidateWorker(mixed_binding, binding)$valid)

    non_finite <- list(
        schema = .OMICS_MEMORY_WORKER_SCHEMA,
        schema_version = .OMICS_MEMORY_WORKER_VERSION,
        status = "passed",
        binding_sha256 = binding$binding_sha256,
        payload_sha256 = binding$payload_sha256,
        truth_sha256 = binding$truth_sha256,
        retained_rss_bytes = Inf
    )
    expect_false(omicsMemoryValidateWorker(non_finite, binding)$valid)

    malformed_adapter <- list(status = "passed")
    expect_error(
        omicsMemoryValidateAdapterRun(malformed_adapter),
        class = "omics_workload_contract_error"
    )

    run_dir <- tempfile("omics-undeclared-file-")
    dir.create(run_dir)
    writeLines("not declared", file.path(run_dir, "rogue.txt"))
    expect_identical(omicsMemoryUndeclaredFiles(run_dir), "rogue.txt")
})
