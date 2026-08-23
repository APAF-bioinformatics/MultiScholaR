.proteomicsWorkloadRepoPath <- function(...) {
    file.path(normalizePath(testthat::test_path("..", ".."), mustWork = TRUE), ...)
}

source(
    .proteomicsWorkloadRepoPath(
        "tools",
        "profiling",
        "omics_workload_contract.R"
    ),
    local = TRUE
)

.proteomicsWorkloadAdapterPath <- .proteomicsWorkloadRepoPath(
    "tools",
    "profiling",
    "omics_workload_proteomics.R"
)

.proteomicsWorkloadContractPaths <- function() {
    paths <- list.files(
        .proteomicsWorkloadRepoPath(
            "tests",
            "testdata",
            "omics-parity",
            "proteomics",
            "workloads"
        ),
        pattern = "[.]json$",
        full.names = TRUE
    )
    sort(paths, method = "radix")
}

.proteomicsWorkloadPrepare <- function(contract, adapter) {
    do.call(RNGkind, as.list(unlist(contract$rng$kind, use.names = FALSE)))
    set.seed(as.integer(contract$rng$seed))
    run_dir <- tempfile("proteomics-workload-test-")
    dir.create(run_dir, recursive = TRUE)
    prepared <- adapter$prepare(list(
        contract = contract,
        run_dir = run_dir,
        repetition = 1L
    ))
    list(run_dir = run_dir, prepared = prepared)
}

test_that("non-DIA proteomics workload contracts bind exact supported tuples", {
    paths <- .proteomicsWorkloadContractPaths()
    expect_length(paths, 3L)

    contracts <- lapply(paths, omicsWorkloadReadContract)
    capability_ids <- vapply(
        contracts,
        \(contract) contract$capability$capability_id,
        character(1)
    )
    expect_setequal(capability_ids, c(
        "proteomics.maxquant.protein.lfq.v1",
        "proteomics.fragpipe.protein.lfq.v1",
        "proteomics.pd_tmt.protein.tmt.v1"
    ))

    for (contract in contracts) {
        expect_identical(contract$privacy$classification, "public_synthetic")
        expect_identical(contract$execution$repetitions, 5L)
        expect_identical(contract$dimensions$feature_count, 10000L)
        expect_identical(contract$dimensions$sample_count, 12L)
        expect_identical(
            contract$adapter$source_sha256,
            omicsWorkloadFileDigest(.proteomicsWorkloadAdapterPath)
        )
        expect_invisible(omicsWorkloadValidateContract(contract))
    }
})

test_that("proteomics adapter reproduces payload, truth, and importer evidence", {
    for (path in .proteomicsWorkloadContractPaths()) {
        contract <- omicsWorkloadReadContract(path)
        adapter <- omicsWorkloadLoadAdapter(
            .proteomicsWorkloadAdapterPath,
            contract
        )
        first <- .proteomicsWorkloadPrepare(contract, adapter)
        second <- .proteomicsWorkloadPrepare(contract, adapter)
        first_prepared <- first$prepared
        second_prepared <- second$prepared

        expect_identical(
            omicsWorkloadFileDigest(first_prepared$payload_path),
            contract$expected_digests$payload_sha256,
            info = basename(path)
        )
        expect_identical(
            omicsWorkloadFileDigest(first_prepared$truth_path),
            contract$expected_digests$truth_sha256,
            info = basename(path)
        )
        expect_identical(
            omicsWorkloadFileDigest(first_prepared$payload_path),
            omicsWorkloadFileDigest(second_prepared$payload_path),
            info = basename(path)
        )
        expect_identical(
            omicsWorkloadFileDigest(first_prepared$truth_path),
            omicsWorkloadFileDigest(second_prepared$truth_path),
            info = basename(path)
        )

        result <- adapter$run(list(
            contract = contract,
            payload_path = first_prepared$payload_path,
            truth_path = first_prepared$truth_path,
            run_dir = first$run_dir,
            package_library = NULL,
            diagnostics_only = FALSE
        ))
        expect_true(result$workflow_evidence$truth_valid, info = basename(path))
        expect_identical(
            result$workflow_evidence$import$rows,
            120000L,
            info = basename(path)
        )
        expect_identical(
            result$query_evidence[[1L]]$rows,
            1200L,
            info = basename(path)
        )
        expect_identical(result$session_evidence$resource_count, 0L)
        expect_identical(result$report_evidence$file_count, 0L)
    }
})

test_that("tuple and seed changes cannot pool proteomics workload evidence", {
    paths <- .proteomicsWorkloadContractPaths()
    first <- omicsWorkloadReadContract(paths[[1L]])
    second <- omicsWorkloadReadContract(paths[[2L]])
    expect_false(identical(omicsWorkloadDigest(first), omicsWorkloadDigest(second)))
    expect_false(identical(
        first$expected_digests$payload_sha256,
        second$expected_digests$payload_sha256
    ))

    changed_seed <- first
    changed_seed$rng$seed <- changed_seed$rng$seed + 1L
    changed_seed$expected_digests <- list(
        payload_sha256 = paste(rep("0", 64L), collapse = ""),
        truth_sha256 = paste(rep("0", 64L), collapse = "")
    )
    adapter <- omicsWorkloadLoadAdapter(
        .proteomicsWorkloadAdapterPath,
        changed_seed
    )
    changed <- .proteomicsWorkloadPrepare(changed_seed, adapter)$prepared
    expect_false(identical(
        omicsWorkloadFileDigest(changed$payload_path),
        first$expected_digests$payload_sha256
    ))
    expect_false(identical(
        omicsWorkloadFileDigest(changed$truth_path),
        first$expected_digests$truth_sha256
    ))
})

test_that("proteomics adapter rejects tampered workload truth", {
    path <- .proteomicsWorkloadContractPaths()[[1L]]
    contract <- omicsWorkloadReadContract(path)
    adapter <- omicsWorkloadLoadAdapter(
        .proteomicsWorkloadAdapterPath,
        contract
    )
    prepared <- .proteomicsWorkloadPrepare(contract, adapter)
    truth <- jsonlite::read_json(
        prepared$prepared$truth_path,
        simplifyVector = TRUE
    )
    truth$long_row_count <- truth$long_row_count + 1L
    jsonlite::write_json(
        truth,
        prepared$prepared$truth_path,
        auto_unbox = TRUE,
        pretty = TRUE,
        digits = 17
    )

    expect_error(
        adapter$run(list(
            contract = contract,
            payload_path = prepared$prepared$payload_path,
            truth_path = prepared$prepared$truth_path,
            run_dir = prepared$run_dir,
            package_library = NULL,
            diagnostics_only = FALSE
        )),
        "does not match declared truth",
        fixed = TRUE
    )
})
