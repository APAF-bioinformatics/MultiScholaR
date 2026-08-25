lipidWorkloadTestRepoPath <- function(...) {
    file.path(normalizePath(testthat::test_path("..", ".."), mustWork = TRUE), ...)
}

lipidWorkloadTestEnvironment <- function() {
    environment <- new.env(parent = .GlobalEnv)
    sys.source(
        lipidWorkloadTestRepoPath(
            "tools", "profiling", "omics_workload_contract.R"
        ),
        envir = environment
    )
    environment
}

lipidWorkloadTestContracts <- function() {
    sort(list.files(
        lipidWorkloadTestRepoPath(
            "tests", "testdata", "omics-parity", "lipidomics", "workloads"
        ),
        pattern = "[.]json$",
        full.names = TRUE
    ))
}

lipidWorkloadTestPrepare <- function(adapter, contract, root) {
    do.call(RNGkind, as.list(unlist(contract$rng$kind, use.names = FALSE)))
    set.seed(as.integer(contract$rng$seed))
    adapter$prepare(list(contract = contract, run_dir = root, repetition = 1L))
}

test_that("lipidomics workload contracts and chemistry are exact", {
    environment <- lipidWorkloadTestEnvironment()
    adapter_path <- lipidWorkloadTestRepoPath(
        "tools", "profiling", "omics_workload_lipidomics.R"
    )
    paths <- lipidWorkloadTestContracts()
    expect_length(paths, 11L)
    ids <- character()
    profiles <- character()
    for (path in paths) {
        contract <- environment$omicsWorkloadReadContract(path)
        adapter <- environment$omicsWorkloadLoadAdapter(adapter_path, contract)
        prepared <- lipidWorkloadTestPrepare(
            adapter,
            contract,
            withr::local_tempdir(pattern = "lipid-workload-")
        )
        expect_identical(
            digest::digest(file = prepared$payload_path, algo = "sha256"),
            contract$expected_digests$payload_sha256
        )
        expect_identical(
            digest::digest(file = prepared$truth_path, algo = "sha256"),
            contract$expected_digests$truth_sha256
        )
        result <- adapter$run(list(
            contract = contract,
            payload_path = prepared$payload_path,
            truth_path = prepared$truth_path,
            run_dir = dirname(prepared$payload_path)
        ))
        expect_identical(result$status, "passed")
        expect_true(result$workflow_evidence$truth_valid)
        if (identical(
            contract$adapter_parameters$source_kind,
            "generated_custom"
        )) {
            expect_identical(contract$capability$input_format, "custom")
            expect_true(prepared$metadata$generated)
            expect_identical(
                sum(!grepl("^LIPID_S[0-9]+$", names(result$retained))),
                9L
            )
            expect_gt(length(result$workflow_evidence$lipid_classes), 0L)
            expect_gt(length(result$workflow_evidence$adducts), 0L)
            expect_gt(length(result$workflow_evidence$ion_modes), 0L)
            expect_true(all(unlist(
                result$workflow_evidence$missing_counts,
                use.names = FALSE
            ) > 0L))
            expect_match(
                result$workflow_evidence$claim_boundary,
                "not vendor or biological validation",
                fixed = TRUE
            )
        } else {
            expect_identical(contract$capability$input_format, "lipidsearch")
            expect_false(prepared$metadata$generated)
            expect_identical(
                result$workflow_evidence$evidence_class,
                "independent_fixture"
            )
            expect_identical(
                result$workflow_evidence$assay,
                contract$adapter_parameters$oracle$assay
            )
            expect_match(
                result$workflow_evidence$claim_boundary,
                "not representative-scale or biological validation",
                fixed = TRUE
            )
        }
        ids <- c(ids, contract$workload_id)
        profiles <- c(profiles, contract$adapter_parameters$profile)
    }
    expect_identical(anyDuplicated(ids), 0L)
    expect_identical(anyDuplicated(profiles), 0L)
})

test_that("lipidomics simulator is deterministic and binding-sensitive", {
    environment <- lipidWorkloadTestEnvironment()
    path <- lipidWorkloadTestContracts()[grepl(
        "mixed-public-ci",
        lipidWorkloadTestContracts()
    )][[1L]]
    contract <- environment$omicsWorkloadReadContract(path)
    adapter <- environment$omicsWorkloadLoadAdapter(
        lipidWorkloadTestRepoPath(
            "tools", "profiling", "omics_workload_lipidomics.R"
        ),
        contract
    )
    first <- lipidWorkloadTestPrepare(
        adapter,
        contract,
        withr::local_tempdir(pattern = "lipid-determinism-a-")
    )
    second <- lipidWorkloadTestPrepare(
        adapter,
        contract,
        withr::local_tempdir(pattern = "lipid-determinism-b-")
    )
    expect_identical(
        digest::digest(file = first$payload_path, algo = "sha256"),
        digest::digest(file = second$payload_path, algo = "sha256")
    )
    changed <- contract
    changed$rng$seed <- changed$rng$seed + 1L
    changed$workload_id <- paste0(changed$workload_id, ".changed")
    third <- lipidWorkloadTestPrepare(
        adapter,
        changed,
        withr::local_tempdir(pattern = "lipid-determinism-c-")
    )
    expect_false(identical(
        digest::digest(file = first$payload_path, algo = "sha256"),
        digest::digest(file = third$payload_path, algo = "sha256")
    ))
    expect_false(identical(
        environment$omicsWorkloadDigest(contract),
        environment$omicsWorkloadDigest(changed)
    ))
})

test_that("private scale mapping is lipidomics-specific and sanitized", {
    environment <- lipidWorkloadTestEnvironment()
    sys.source(
        lipidWorkloadTestRepoPath(
            "tools", "profiling", "omics_workload_private_scale.R"
        ),
        envir = environment
    )
    private_path <- lipidWorkloadTestRepoPath("data", "cotton_report.tsv")
    testthat::skip_if_not(file.exists(private_path), "local private input absent")
    before <- file.info(private_path)[c("size", "mtime")]
    manifest <- environment$inspectPrivateScaleTsv(
        private_path,
        paste(rep("local-lipid-scale-salt", 2L), collapse = ":")
    )
    after <- file.info(private_path)[c("size", "mtime")]
    expect_identical(before, after)
    expect_identical(names(manifest), c(
        "row_count", "column_count", "byte_size", "salted_source_fingerprint"
    ))
    mapping <- environment$mapPrivateScaleToLipidomics(manifest)
    expect_identical(mapping$sample_count, 12L)
    expect_identical(length(mapping$assay_mix), 3L)
    expect_identical(mapping$source_report_column_count, 60L)
    expect_false(mapping$sample_count == manifest$column_count)
    expect_match(mapping$mapping_policy, "bounded_lipids", fixed = TRUE)
    expect_false(any(grepl(
        "path|fasta|sequence|accession|header",
        c(names(manifest), names(mapping)),
        ignore.case = TRUE
    )))
})

test_that("private lipidomics scale admission fails closed", {
    environment <- lipidWorkloadTestEnvironment()
    sys.source(
        lipidWorkloadTestRepoPath(
            "tools", "profiling", "omics_workload_private_scale.R"
        ),
        envir = environment
    )
    path <- withr::local_tempfile(fileext = ".tsv")
    writeLines(c("a\tb", "1\t2", "3"), path)
    salt <- paste(rep("lipid-private-scale-salt", 2L), collapse = ":")
    expect_error(
        environment$inspectPrivateScaleTsv(path, salt),
        class = "omics_private_scale_nonrectangular"
    )
    writeLines(c("a\tb", "1\t2"), path)
    expect_error(
        environment$inspectPrivateScaleTsv(path, ""),
        class = "omics_private_scale_unsalted"
    )
    expect_error(
        environment$inspectPrivateScaleTsv(paste0(path, ".missing"), salt),
        class = "omics_private_scale_unreadable"
    )

    original_state <- environment$privateScaleFileState
    state_calls <- 0L
    environment$privateScaleFileState <- function(candidate) {
        state_calls <<- state_calls + 1L
        state <- original_state(candidate)
        if (state_calls > 1L) state$modified <- state$modified + 1
        state
    }
    expect_error(
        environment$inspectPrivateScaleTsv(path, salt),
        class = "omics_private_scale_changed"
    )
    expect_error(
        environment$mapPrivateScaleToLipidomics(list(row_count = 100L)),
        class = "omics_private_scale_unsafe_manifest"
    )
})

test_that("private-calibrated evidence is separate and path-free", {
    environment <- lipidWorkloadTestEnvironment()
    private_contract <- "/tmp/multischolar-lipid-private/private-contract.json"
    private_evidence <-
        "/tmp/multischolar-lipid-private-baseline/baseline-result.json"
    testthat::skip_if_not(
        file.exists(private_contract) && file.exists(private_evidence),
        "opt-in private lipid calibration evidence absent"
    )
    private <- environment$omicsWorkloadReadContract(private_contract)
    public <- environment$omicsWorkloadReadContract(
        lipidWorkloadTestContracts()[grepl(
            "mixed-public-representative",
            lipidWorkloadTestContracts()
        )][[1L]]
    )
    expect_identical(private$privacy$classification, "private_scale_only")
    expect_false(identical(private$workload_id, public$workload_id))
    private_science <- private$adapter_parameters
    public_science <- public$adapter_parameters
    scale_fields <- c(
        "profile", "scale_mapping_policy", "source_report_column_count",
        "source_byte_class"
    )
    private_science[scale_fields] <- NULL
    public_science[scale_fields] <- NULL
    expect_identical(private_science, public_science)
    expect_identical(private$dimensions$sample_count, 12L)
    expect_identical(private$dimensions$assay_count, 3L)
    expect_false(private$dimensions$sample_count ==
        private$privacy$scale_metadata$column_count)
    evidence <- jsonlite::read_json(private_evidence, simplifyVector = FALSE)
    expect_identical(evidence$status, "passed")
    expect_identical(evidence$summary$completed, 5L)
    retained <- paste(readLines(private_evidence, warn = FALSE), collapse = "\n")
    local_source <- lipidWorkloadTestRepoPath("data", "cotton_report.tsv")
    expect_false(grepl(basename(local_source), retained, fixed = TRUE))
    expect_false(grepl(
        normalizePath(dirname(local_source), mustWork = TRUE),
        retained,
        fixed = TRUE
    ))
    expect_false(grepl("[.](fasta|fa|crdownload)", retained, ignore.case = TRUE))
})
