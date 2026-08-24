metabWorkloadTestRepoPath <- function(...) {
    file.path(normalizePath(testthat::test_path("..", ".."), mustWork = TRUE), ...)
}

metabWorkloadTestSourceHarness <- function() {
    source(metabWorkloadTestRepoPath(
        "tools", "profiling", "omics_workload_contract.R"
    ), local = .GlobalEnv)
    source(metabWorkloadTestRepoPath(
        "tools", "profiling", "omics_workload_private_scale.R"
    ), local = .GlobalEnv)
}

metabWorkloadTestContracts <- function() {
    sort(list.files(
        metabWorkloadTestRepoPath(
            "tests", "testdata", "omics-parity", "metabolomics", "workloads"
        ),
        pattern = "[.]json$",
        full.names = TRUE
    ))
}

metabWorkloadTestPrepare <- function(adapter, contract, root) {
    do.call(RNGkind, as.list(unlist(contract$rng$kind, use.names = FALSE)))
    set.seed(as.integer(contract$rng$seed))
    adapter$prepare(list(contract = contract, run_dir = root, repetition = 1L))
}

test_that("metabolomics workload contracts and adapter are exact", {
    metabWorkloadTestSourceHarness()
    adapter_path <- metabWorkloadTestRepoPath(
        "tools", "profiling", "omics_workload_metabolomics.R"
    )
    contracts <- metabWorkloadTestContracts()
    expect_length(contracts, 5L)
    workload_ids <- character()
    profiles <- character()
    for (path in contracts) {
        contract <- omicsWorkloadReadContract(path)
        adapter <- omicsWorkloadLoadAdapter(adapter_path, contract)
        expect_identical(adapter$adapter_id, "multischolar.metabolomics.custom.v1")
        expect_identical(contract$capability$input_format, "custom")
        expect_identical(contract$capability$omic_type, "metabolomics")
        expect_identical(
            sort(names(contract$missingness$mechanisms)),
            sort(c("MCAR", "MAR", "MNAR", "left_censored"))
        )
        expect_gt(contract$adapter_parameters$internal_standard_fraction, 0)
        expect_gt(contract$adapter_parameters$heteroscedasticity, 0)
        expect_gt(contract$adapter_parameters$batch_log2, 0)
        root <- withr::local_tempdir(pattern = "metab-workload-")
        prepared <- metabWorkloadTestPrepare(adapter, contract, root)
        expect_identical(
            digest::digest(
                file = prepared$payload_path,
                algo = "sha256",
                serialize = FALSE
            ),
            contract$expected_digests$payload_sha256
        )
        expect_identical(
            digest::digest(
                file = prepared$truth_path,
                algo = "sha256",
                serialize = FALSE
            ),
            contract$expected_digests$truth_sha256
        )
        result <- adapter$run(list(
            contract = contract,
            payload_path = prepared$payload_path,
            truth_path = prepared$truth_path,
            run_dir = root
        ))
        expect_identical(result$status, "passed")
        expect_true(result$workflow_evidence$truth_valid)
        expect_match(
            result$workflow_evidence$claim_boundary,
            "not vendor or biological validation",
            fixed = TRUE
        )
        expect_gt(result$workflow_evidence$internal_standard_count, 0L)
        expect_gt(result$workflow_evidence$duplicate_id_count, 0L)
        expect_true(all(unlist(
            result$workflow_evidence$missing_counts,
            use.names = FALSE
        ) > 0L))
        expect_length(result$query_evidence, 1L)
        expect_lte(
            result$query_evidence[[1L]]$rows,
            contract$adapter_parameters$bounded_query_features * 2L
        )
        workload_ids <- c(workload_ids, contract$workload_id)
        profiles <- c(profiles, contract$adapter_parameters$profile)
    }
    expect_identical(anyDuplicated(workload_ids), 0L)
    expect_identical(anyDuplicated(profiles), 0L)
})

test_that("metabolomics simulator is deterministic and binding-sensitive", {
    metabWorkloadTestSourceHarness()
    path <- metabWorkloadTestContracts()[grepl(
        "mixed-public-ci",
        metabWorkloadTestContracts()
    )][[1L]]
    contract <- omicsWorkloadReadContract(path)
    adapter <- omicsWorkloadLoadAdapter(
        metabWorkloadTestRepoPath(
            "tools", "profiling", "omics_workload_metabolomics.R"
        ),
        contract
    )
    first <- metabWorkloadTestPrepare(
        adapter,
        contract,
        withr::local_tempdir(pattern = "metab-determinism-a-")
    )
    second <- metabWorkloadTestPrepare(
        adapter,
        contract,
        withr::local_tempdir(pattern = "metab-determinism-b-")
    )
    expect_identical(
        digest::digest(file = first$payload_path, algo = "sha256"),
        digest::digest(file = second$payload_path, algo = "sha256")
    )
    expect_identical(
        digest::digest(file = first$truth_path, algo = "sha256"),
        digest::digest(file = second$truth_path, algo = "sha256")
    )
    changed <- contract
    changed$rng$seed <- changed$rng$seed + 1L
    changed$workload_id <- paste0(changed$workload_id, ".changed-seed")
    changed_prepared <- metabWorkloadTestPrepare(
        adapter,
        changed,
        withr::local_tempdir(pattern = "metab-determinism-c-")
    )
    expect_false(identical(
        digest::digest(file = first$payload_path, algo = "sha256"),
        digest::digest(file = changed_prepared$payload_path, algo = "sha256")
    ))
    expect_false(identical(
        omicsWorkloadDigest(contract),
        omicsWorkloadDigest(changed)
    ))
})

test_that("private DIA scale admission is sanitized and metabolomics-specific", {
    metabWorkloadTestSourceHarness()
    private_path <- metabWorkloadTestRepoPath("data", "cotton_report.tsv")
    testthat::skip_if_not(file.exists(private_path), "local private scale input absent")
    salt <- paste(rep("local-closeout-salt", 2L), collapse = ":")
    before <- file.info(private_path)[c("size", "mtime")]
    manifest <- inspectPrivateScaleTsv(private_path, salt)
    after <- file.info(private_path)[c("size", "mtime")]
    expect_identical(before, after)
    expect_identical(names(manifest), c(
        "row_count", "column_count", "byte_size",
        "salted_source_fingerprint"
    ))
    expect_identical(manifest$row_count, 80170L)
    expect_identical(manifest$column_count, 60L)
    expect_identical(manifest$byte_size, 71000019)
    expect_match(manifest$salted_source_fingerprint, "^[a-f0-9]{64}$")
    expect_false(identical(
        manifest$salted_source_fingerprint,
        digest::digest(file = private_path, algo = "sha256")
    ))
    mapping <- mapPrivateScaleToMetabolomics(manifest)
    expect_identical(mapping$sample_count, 12L)
    expect_identical(length(mapping$assay_mix), 3L)
    expect_identical(mapping$source_report_column_count, 60L)
    expect_false(mapping$feature_count == manifest$row_count)
    expect_false(mapping$sample_count == manifest$column_count)
    expect_identical(
        mapping$mapping_policy,
        "rows_to_bounded_features_columns_io_only_v1"
    )
    fasta_path <- metabWorkloadTestRepoPath(
        "data", "Unconfirmed 265895.crdownload"
    )
    expect_false(any(grepl(
        basename(fasta_path),
        c(names(manifest), names(mapping)),
        fixed = TRUE
    )))
})

test_that("private scale admission rejects unsafe inputs without fallback", {
    metabWorkloadTestSourceHarness()
    path <- withr::local_tempfile(fileext = ".tsv")
    writeLines(c("a\tb", "1\t2", "3"), path)
    expect_error(
        inspectPrivateScaleTsv(path, paste(rep("salt", 5L), collapse = ":")),
        class = "omics_private_scale_nonrectangular"
    )
    expect_error(
        inspectPrivateScaleTsv(path, ""),
        class = "omics_private_scale_unsalted"
    )
    expect_error(
        inspectPrivateScaleTsv(paste0(path, ".missing"), paste(rep("salt", 5L), collapse = ":")),
        class = "omics_private_scale_unreadable"
    )
    unsafe <- list(
        row_count = 10L,
        column_count = 4L,
        byte_size = 100,
        salted_source_fingerprint = strrep("a", 64L),
        source_path = "/private/source.tsv"
    )
    expect_error(
        mapPrivateScaleToMetabolomics(unsafe),
        class = "omics_private_scale_unsafe_manifest"
    )
})

test_that("private-calibrated evidence differs only in declared scale", {
    metabWorkloadTestSourceHarness()
    private_contract_path <-
        "/tmp/multischolar-metab-private/private-contract.json"
    private_evidence_path <-
        "/tmp/multischolar-metab-private-baseline/baseline-result.json"
    testthat::skip_if_not(
        file.exists(private_contract_path) && file.exists(private_evidence_path),
        "opt-in private calibration evidence absent"
    )
    private <- omicsWorkloadReadContract(private_contract_path)
    public <- omicsWorkloadReadContract(metabWorkloadTestContracts()[grepl(
        "mixed-public-representative",
        metabWorkloadTestContracts()
    )][[1L]])
    expect_identical(private$privacy$classification, "private_scale_only")
    expect_false(identical(private$workload_id, public$workload_id))
    expect_false(identical(
        omicsWorkloadDigest(private),
        omicsWorkloadDigest(public)
    ))
    private_science <- private$adapter_parameters
    public_science <- public$adapter_parameters
    private_science$profile <- NULL
    public_science$profile <- NULL
    expect_identical(private_science, public_science)
    expect_identical(private$dimensions$sample_count, 12L)
    expect_identical(private$dimensions$assay_count, 3L)
    expect_false(private$dimensions$sample_count ==
        private$privacy$scale_metadata$column_count)
    expect_silent(omicsWorkloadForbiddenPrivateFields(private))
    evidence <- jsonlite::read_json(
        private_evidence_path,
        simplifyVector = FALSE
    )
    expect_identical(evidence$status, "passed")
    expect_identical(evidence$summary$completed, 5L)
    expect_identical(
        evidence$runs[[1L]]$binding$privacy_classification,
        "private_scale_only"
    )
    retained <- paste(readLines(private_evidence_path, warn = FALSE), collapse = "\n")
    expect_false(grepl("cotton_report", retained, fixed = TRUE))
    expect_false(grepl("Unconfirmed 265895", retained, fixed = TRUE))
    expect_false(grepl("/home/doktersmol", retained, fixed = TRUE))
})
