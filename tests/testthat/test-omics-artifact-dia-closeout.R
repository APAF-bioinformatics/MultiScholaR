diaArtifact020SkipDependencies <- function() {
    for (package in c("arrow", "DBI", "duckdb", "filelock", "processx", "ps")) {
        testthat::skip_if_not_installed(package)
    }
}

diaArtifact020Runner <- function() {
    omicsParityRepoPath(
        "tools",
        "profiling",
        "run_omics_artifact_baseline.R"
    )
}

diaArtifact020LoadHelpers <- function(envir = parent.frame()) {
    source(
        omicsParityRepoPath(
            "tools",
            "profiling",
            "omics_artifact_closeout_helpers.R"
        ),
        local = envir
    )
}

diaArtifact020RunCloseout <- function(
    scenario,
    output_dir,
    include_private = FALSE,
    env = character()
) {
    processx::run(
        file.path(R.home("bin"), "Rscript"),
        c(
            "--vanilla",
            diaArtifact020Runner(),
            "--scenario",
            scenario,
            "--backend",
            "paired",
            "--repetitions",
            "1",
            "--diagnostics",
            "false",
            "--include-private",
            tolower(as.character(include_private)),
            "--output-dir",
            output_dir
        ),
        wd = tempdir(),
        env = env,
        timeout = 180000,
        error_on_status = FALSE,
        echo = FALSE
    )
}

diaArtifact020MockRun <- function(backend, repetition) {
    artifact <- identical(backend, "artifact")
    query_stage <- list(
        stage_id = "bounded_query",
        results = list(list(output_sha256 = "query-digest"))
    )
    list(
        status = "passed",
        scenario_id = "representative",
        backend = backend,
        repetition = as.integer(repetition),
        metrics = list(
            peak_tree_rss_bytes = if (artifact) 60 else 100,
            retained_tree_rss_bytes = if (artifact) 50 else 100,
            elapsed_seconds = if (artifact) 12 else 10,
            committed_input_bytes = 100,
            committed_artifact_bytes = if (artifact) 120 else 0,
            peak_artifact_disk_bytes = if (artifact) 200 else 0,
            final_file_count = if (artifact) 6 else 2,
            bounded_query_p95_seconds = if (artifact) 0.3 else 0.2
        ),
        worker = list(
            status = "passed",
            observed_summary = list(output_sha256 = "summary-digest"),
            stages = if (artifact) {
                list(
                    list(stage_id = "import"),
                    list(stage_id = "artifact_persist"),
                    query_stage
                )
            } else {
                list(list(stage_id = "import"), query_stage)
            }
        )
    )
}

test_that("DIA E2E manifest freezes the exact artifact closeout scope", {
    manifest <- jsonlite::read_json(
        omicsParityRepoPath("tests", "testdata", "e2e", "manifest.json"),
        simplifyVector = FALSE
    )
    lanes <- manifest$lanes
    names(lanes) <- vapply(lanes, `[[`, character(1), "lane_id")
    closeout <- lanes$prot_dia$artifact_closeout

    expect_identical(
        closeout$capability_id,
        "proteomics.diann.peptide.dia.v1"
    )
    expect_identical(
        closeout$performance_scenario,
        "dia_private_representative"
    )
    expect_identical(closeout$minimum_paired_repetitions, 5L)
    expect_false(closeout$generated_evidence_authorizes_promotion)
    expect_null(lanes$prot_dia_limpa$artifact_closeout)
})

test_that("DIA closeout promotion evaluation is paired and fail closed", {
    diaArtifact020LoadHelpers()
    gates <- omicsParityReadManifest()$release_gates
    runs <- unlist(lapply(seq_len(5L), function(repetition) {
        list(
            diaArtifact020MockRun("memory", repetition),
            diaArtifact020MockRun("artifact", repetition)
        )
    }), recursive = FALSE)
    scenario <- list(scenario_id = "representative", kind = "optional_private")

    passed <- baselinePromotionEvaluation(runs, list(scenario), gates)
    expect_true(passed$authorized)
    expect_false(passed$generated_evidence_authorizes_promotion)
    expect_true(all(vapply(
        passed$scenarios[[1L]]$gates,
        `[[`,
        logical(1),
        "passed"
    )))

    insufficient <- baselinePromotionEvaluation(runs[seq_len(8L)], list(scenario), gates)
    expect_false(insufficient$authorized)
    expect_identical(
        insufficient$scenarios[[1L]]$reason,
        "insufficient_paired_repetitions"
    )

    generated <- scenario
    generated$kind <- "generated_scaling"
    generated_result <- baselinePromotionEvaluation(runs, list(generated), gates)
    expect_false(generated_result$authorized)
    expect_identical(
        generated_result$scenarios[[1L]]$reason,
        "non_representative_evidence"
    )

    incomplete <- runs
    incomplete[[2L]]$metrics$bounded_query_p95_seconds <- NA_real_
    incomplete_result <- baselinePromotionEvaluation(
        incomplete,
        list(scenario),
        gates
    )
    expect_false(incomplete_result$authorized)
    query_gate <- Filter(function(gate) {
        identical(gate$gate_id, "maximum_bounded_query_p95_seconds")
    }, incomplete_result$scenarios[[1L]]$gates)[[1L]]
    expect_identical(query_gate$valid_pairs, 4L)
    expect_false(query_gate$passed)
})

test_that("paired DIA closeout preserves exact import and bounded query results", {
    skip_on_cran()
    diaArtifact020SkipDependencies()
    output_dir <- tempfile("dia-artifact-020-paired-")
    process <- diaArtifact020RunCloseout("dia_canonical_import", output_dir)
    expect_identical(process$status, 0L, info = process$stderr)

    diaArtifact020LoadHelpers()
    result <- baselineValidateCloseoutResult(
        file.path(output_dir, "baseline-result.json")
    )
    expect_setequal(vapply(result$runs, `[[`, character(1), "backend"), c(
        "memory",
        "artifact"
    ))
    expect_true(all(vapply(result$runs, function(run) {
        isTRUE(run$worker$oracle_comparison$equal)
    }, logical(1))))
    evaluation <- result$promotion_evaluation$scenarios[[1L]]
    expect_true(evaluation$scientific_parity$passed)
    expect_false(evaluation$promotion_eligible)
    expect_false(result$promotion_evaluation$authorized)

    artifact <- Filter(function(run) {
        identical(run$backend, "artifact")
    }, result$runs)[[1L]]
    expect_false(artifact$worker$source_table_retained)
    expect_identical(artifact$worker$resource_evidence$staging_entries, 0L)
    expect_identical(
        artifact$worker$resource_evidence$duckdb_temporary_entries,
        0L
    )
    descriptor <- artifactDiaWorkflowDescriptor()
    expect_identical(descriptor$certification$status, "dual_write")
    expect_false(descriptor$certification$auto_eligible)
    expect_true(all(vapply(descriptor$stages, function(stage) {
        identical(stage$maximum_rollout, "dual_write")
    }, logical(1))))
})

test_that("private paired evidence removes payloads and biological identifiers", {
    skip_on_cran()
    diaArtifact020SkipDependencies()
    fixture <- omicsParityRepoPath(
        "tests",
        "testdata",
        "e2e",
        "prot_dia",
        "report.tsv"
    )
    private_dir <- tempfile("private-dia-closeout-source-")
    dir.create(private_dir)
    private_path <- file.path(private_dir, "private-dia-closeout-marker.tsv")
    expect_true(file.copy(fixture, private_path))
    withr::defer(unlink(private_dir, recursive = TRUE, force = TRUE))
    output_dir <- tempfile("dia-artifact-020-private-")
    process <- diaArtifact020RunCloseout(
        "dia_private_representative",
        output_dir,
        include_private = TRUE,
        env = c(
            MULTISCHOLAR_DIA_BASELINE_FIXTURE = private_path,
            MULTISCHOLAR_BASELINE_FINGERPRINT_SALT = "test-only-closeout-salt"
        )
    )
    expect_identical(process$status, 0L, info = process$stderr)

    diaArtifact020LoadHelpers()
    result <- baselineValidateCloseoutResult(
        file.path(output_dir, "baseline-result.json")
    )
    evaluation <- result$promotion_evaluation$scenarios[[1L]]
    expect_true(evaluation$promotion_eligible)
    expect_identical(evaluation$paired_repetitions, 1L)
    expect_identical(evaluation$reason, "insufficient_paired_repetitions")
    expect_false(result$promotion_evaluation$authorized)

    artifact <- Filter(function(run) {
        identical(run$backend, "artifact")
    }, result$runs)[[1L]]
    expect_true(artifact$worker$private_artifact_payload_retained_at_worker_exit)
    expect_false(artifact$worker$private_artifact_payload_retained_in_evidence)
    expect_false(dir.exists(file.path(
        output_dir,
        "runs",
        "dia_private_representative-artifact-r01",
        "committed",
        "project"
    )))

    persisted <- list.files(output_dir, recursive = TRUE, full.names = TRUE)
    persisted_text <- unlist(lapply(persisted, function(path) {
        paste(readLines(path, warn = FALSE), collapse = "\n")
    }), use.names = FALSE)
    expect_false(any(grepl(private_path, persisted_text, fixed = TRUE)))
    expect_false(any(grepl("WT_1|P00001|PEP1AK", persisted_text)))
})
