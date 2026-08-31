diaRepairTestMeasurement <- function(metrics, proof) {
    list(
        status = "passed",
        publication_certifiable = TRUE,
        cleanup = list(valid = TRUE),
        metrics = c(list(peak_swap_bytes = 0), metrics),
        worker = list(
            scientific_proof = proof,
            complete_payload_returned = FALSE,
            payload_free_state_manager = TRUE
        )
    )
}

diaRepairTestRecords <- function(elapsed_ratio = 0.9) {
    proof <- list(
        salted_state_digest = strrep("a", 64L),
        s4_class = "PeptideQuantitativeData",
        valid_s4 = TRUE,
        source_fields_released = TRUE,
        current_state = "raw_data_s4"
    )
    pre <- list(
        peak_charged_memory_bytes = 100,
        retained_charged_memory_bytes = 100,
        elapsed_seconds = 100,
        primary_work_units_per_wall_second = 100,
        primary_work_units_per_cpu_second = 100,
        phase_cpu_seconds = 100,
        peak_logical_disk_bytes = 100,
        final_logical_disk_bytes = 100
    )
    candidate <- list(
        peak_charged_memory_bytes = 80,
        retained_charged_memory_bytes = 80,
        elapsed_seconds = 100 * elapsed_ratio,
        primary_work_units_per_wall_second = 110,
        primary_work_units_per_cpu_second = 110,
        phase_cpu_seconds = 90,
        peak_logical_disk_bytes = 100,
        final_logical_disk_bytes = 100
    )
    records <- list()
    for (index in seq_len(36L)) {
        pair_id <- sprintf("dia-repair-pair-%03d", index)
        records[[length(records) + 1L]] <- list(
            pair_id = pair_id,
            arm = "pre_repair_artifact",
            measurement = diaRepairTestMeasurement(pre, proof)
        )
        records[[length(records) + 1L]] <- list(
            pair_id = pair_id,
            arm = "candidate_artifact",
            measurement = diaRepairTestMeasurement(candidate, proof)
        )
    }
    records
}

test_that("DIA repair successor preserves gates and corrects owned load", {
    gates <- diaRepairReadGates()

    expect_identical(gates$status, "frozen_after_validated_owned_load_floor")
    expect_identical(
        gates$comparison$predecessor_gate_id,
        "multischolar.dia_commit_repair.2026-08-31.v3"
    )
    expect_identical(
        gates$comparison$predecessor_gate_sha256,
        publicationFileDigest(gates$comparison$predecessor_gate_path)
    )
    expect_identical(
        gates$comparison$phase_start,
        "import_artifact_staging_start"
    )
    expect_identical(
        gates$comparison$peak_scope,
        "complete_owned_cgroup_lifecycle"
    )
    expect_identical(gates$design$required_pairs, 36L)
    expect_identical(gates$design$required_sessions, 3L)
    expect_identical(
        gates$design$host_safety$maximum_prelaunch_thermal_celsius,
        70L
    )
    expect_identical(gates$design$maximum_idle_cpu_activity_seconds, 0.001)
    expect_identical(
        gates$design$host_safety$owned_workload_load_allowance,
        4L
    )
    expect_identical(gates$comparison$manual_garbage_collection_allowed, FALSE)
    expect_false(gates$decision$automatic_policy_authority)
    expect_false(gates$decision$publication_authority)
    expect_setequal(
        vapply(gates$gates, `[[`, character(1), "metric"),
        c(
            "peak_charged_memory_bytes",
            "retained_charged_memory_bytes",
            "elapsed_seconds",
            "primary_work_units_per_wall_second",
            "primary_work_units_per_cpu_second",
            "phase_cpu_seconds",
            "peak_logical_disk_bytes",
            "final_logical_disk_bytes"
        )
    )
})

test_that("DIA repair resume accepts only the bound equivalent predecessor", {
    gates <- diaRepairReadGates()
    predecessor <- list(
        gate_id = gates$comparison$predecessor_gate_id,
        sha256 = gates$comparison$predecessor_gate_sha256
    )

    expect_true(diaRepairPredecessorGateValid(gates))
    expect_true(diaRepairResumeGateBindingValid(predecessor, gates))
    predecessor$sha256 <- strrep("0", 64L)
    expect_false(diaRepairResumeGateBindingValid(predecessor, gates))
})

test_that("DIA repair readiness reserves absolute runtime load headroom", {
    gates <- diaRepairReadGates()
    host <- list(cpu = list(logical_cores = 16L))

    expect_identical(diaRepairPrelaunchMaximumLoad(host, gates), 4)
})

test_that("DIA repair schedule is pair session and order balanced", {
    schedule <- diaRepairSchedule()

    expect_identical(nrow(schedule), 72L)
    expect_setequal(unique(schedule$session_id), sprintf("session-%02d", 1:3))
    expect_true(all(table(schedule$pair_id) == 2L))
    expect_true(all(table(schedule$pair_id, schedule$arm) == 1L))
    expect_true(all(table(schedule$session_id) == 24L))
    first <- schedule[schedule$sequence_in_pair == 1L, ]
    expect_true(all(table(first$session_id, first$arm) == 6L))
})

test_that("DIA repair scientific digests ignore representation details", {
    value <- data.frame(
        id = sprintf("P%04d", seq_len(100L)),
        abundance = seq_len(100L) / 7,
        stringsAsFactors = FALSE
    )
    reordered <- value
    attributes(reordered) <- attributes(reordered)[c(
        "class",
        "names",
        "row.names"
    )]

    expect_identical(value, reordered)
    expect_identical(
        diaRepairStableDigest(value),
        diaRepairStableDigest(reordered)
    )
})

test_that("DIA repair evaluation requires every frozen gate", {
    passing <- diaRepairEvaluate(diaRepairTestRecords())
    failing <- diaRepairEvaluate(diaRepairTestRecords(elapsed_ratio = 1.2))

    expect_identical(passing$status, "passed")
    expect_true(passing$may_start_omics_art_071)
    expect_true(all(vapply(
        passing$numerical_gates,
        `[[`,
        logical(1),
        "passed"
    )))
    expect_identical(failing$status, "failed")
    expect_false(failing$may_start_omics_art_071)
    elapsed <- Filter(function(gate) {
        identical(gate$gate_id, "elapsed_time")
    }, failing$numerical_gates)[[1L]]
    expect_false(elapsed$passed)
})

test_that("DIA repair evaluation excludes every incomplete pair endpoint", {
    records <- diaRepairTestRecords()
    records[[72L]]$measurement$status <- "failed"
    records[[72L]]$measurement$publication_certifiable <- FALSE

    evaluation <- diaRepairEvaluate(records)
    endpoint_pairs <- vapply(evaluation$numerical_gates, function(gate) {
        gate$summary$pairs
    }, integer(1))

    expect_identical(evaluation$valid_pairs, 35L)
    expect_true(all(endpoint_pairs == 35L))
    expect_false(any(vapply(
        evaluation$numerical_gates,
        `[[`,
        logical(1),
        "passed"
    )))
})

test_that("DIA repair resume retains only complete valid pairs", {
    records <- diaRepairTestRecords()[seq_len(68L)]
    records[[68L]]$measurement$status <- "failed"
    records[[68L]]$measurement$publication_certifiable <- FALSE

    retained <- diaRepairCompletePairRecords(records)
    pair_ids <- unique(vapply(retained, `[[`, character(1), "pair_id"))

    expect_length(retained, 66L)
    expect_identical(pair_ids, sprintf("dia-repair-pair-%03d", 1:33))
})

test_that("DIA repair evidence code cannot hide memory with manual GC", {
    paths <- testthat::test_path(
        "..",
        "..",
        "tools",
        "profiling",
        c(
            "omics_dia_commit_repair.R",
            "omics_dia_commit_repair_proof.R",
            "run_omics_dia_commit_repair.R"
        )
    )
    sources <- unlist(lapply(paths, readLines, warn = FALSE), use.names = FALSE)

    expect_false(any(grepl("gc\\s*\\(", sources, perl = TRUE)))
    expect_true(all(lengths(lapply(paths, readLines, warn = FALSE)) < 1000L))
    expect_true(all(nchar(sources, type = "width") <= 100L))
})
