publicationMeasureSpecFixture <- function() {
    protocol_path <- "tests/testdata/omics-performance/protocol-v1.json"
    protocol <- publicationReadJson(protocol_path)
    cache <- publicationReadJson(publicationCacheContractPath())
    digest <- paste(rep("a", 64L), collapse = "")
    spec <- list(
        schema = "multischolar.omics_publication_measure_spec",
        schema_version = "1.0.0",
        measure_spec_id = "fixture-measure-spec",
        protocol_binding = list(
            protocol_id = protocol$protocol_id,
            path = protocol_path,
            sha256 = publicationFileDigest(protocol_path)
        ),
        metric_contract = publicationMetricContractBinding(),
        estimand_binding = list(
            estimand_id = "import_and_settle",
            estimands_sha256 = digest
        ),
        schedule_binding = list(
            schedule_id = "schedule-1",
            schedule_sha256 = digest,
            slot_id = "slot-1",
            slot_sha256 = digest,
            work_binding_sha256 = NULL
        ),
        source_binding = list(
            project_id = "project-1",
            source_id = "source-1",
            source_sha256 = digest,
            input_sha256 = digest,
            exact_input_bytes = 1000,
            evidence_class = "public_real",
            private_source = FALSE,
            private_values_retained = FALSE,
            privacy_receipt_sha256 = digest
        ),
        candidate_binding = list(revision = paste(rep("b", 40L), collapse = ""),
            package_sha256 = digest),
        cache_evidence = list(
            cache_contract_id = cache$cache_contract_id,
            stratum = "standardized_warm_input",
            input_sha256 = digest,
            pre_read_complete = TRUE,
            page_cache_drop_attempted = FALSE,
            page_cache_drop_verified = FALSE,
            page_cache_drop_method = NULL,
            page_cache_drop_receipt_sha256 = NULL,
            primary_comparison_eligible = TRUE
        ),
        work = list(
            primary_work_unit_id = "validated_input_byte",
            expected_primary_work_count = 1000,
            exact_input_bytes = 1000
        ),
        command = file.path(R.home("bin"), "Rscript"),
        arguments = list("--vanilla", "worker.R"),
        run_dir = "/tmp/fixture-run",
        execution = list(
            sampling_interval_ms = 20,
            disk_sampling_interval_ms = 500,
            timeout_seconds = 60,
            retained_dwell_seconds = 5,
            retained_window_seconds = 2,
            retention_acknowledgement = "fifo_v1",
            maximum_boundary_bracket_seconds = 0.5,
            retained_boundary_tolerance_ms = 100
        ),
        environment = list(
            OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
            MKL_NUM_THREADS = "1", ARROW_NUM_THREADS = "1",
            DUCKDB_THREADS = "1", TZ = "UTC"
        ),
        host_preflight = list(
            schema = "multischolar.omics_publication_host_preflight",
            certified = TRUE
        ),
        host_preflight_sha256 = NULL,
        safety_limits = list(
            minimum_available_memory_bytes = 1,
            minimum_available_disk_bytes = 1,
            maximum_load = 1,
            maximum_thermal_celsius = 85,
            baseline_cpu_frequency_hz = 1,
            maximum_cpu_frequency_drift_fraction = 0.05,
            maximum_run_allocated_disk_bytes = 1
        )
    )
    spec$host_preflight_sha256 <- publicationObjectDigest(
        spec$host_preflight
    )
    spec$schedule_binding$work_binding_sha256 <- publicationObjectDigest(
        spec$work
    )
    spec
}

test_that("publication measure specs bind protocol metrics source and work", {
    spec <- publicationMeasureSpecFixture()
    expect_silent(publicationValidateMeasureSpec(spec))

    drift <- publicationGovernanceCopy(spec)
    drift$source_binding$source_sha256 <- paste(rep("0", 63L), collapse = "")
    expect_error(
        publicationValidateMeasureSpec(drift),
        class = "multischolar_publication_binding_error"
    )

    threaded <- publicationGovernanceCopy(spec)
    threaded$environment$DUCKDB_THREADS <- "4"
    expect_error(
        publicationValidateMeasureSpec(threaded),
        class = "multischolar_publication_measure_spec_error"
    )

    input_drift <- publicationGovernanceCopy(spec)
    input_drift$source_binding$exact_input_bytes <- 999
    expect_error(
        publicationValidateMeasureSpec(input_drift),
        class = "multischolar_publication_binding_error"
    )

    privacy_drift <- publicationGovernanceCopy(spec)
    privacy_drift$source_binding$private_values_retained <- TRUE
    expect_error(
        publicationValidateMeasureSpec(privacy_drift),
        class = "multischolar_publication_binding_error"
    )
})

test_that("measure specs accept dynamic observed safety limits additively", {
    spec <- publicationMeasureSpecFixture()
    spec$host_preflight$schema_version <- "1.1.0"
    spec$host_preflight$power_policy <- list(
        mode = "observe_do_not_modify",
        governors = list("powersave"),
        mutation_allowed = FALSE
    )
    spec$host_preflight_sha256 <- publicationObjectDigest(
        spec$host_preflight
    )
    spec$safety_limits <- list(
        safety_mode = "dynamic_observed_v1",
        minimum_available_memory_bytes = 1,
        minimum_available_disk_bytes = 1,
        maximum_load = 1,
        maximum_thermal_celsius = 85,
        baseline_governors = list("powersave"),
        frequency_telemetry_required = TRUE,
        governor_stability_required = TRUE,
        maximum_run_allocated_disk_bytes = 1
    )
    expect_silent(publicationValidateMeasureSpec(spec))

    changed <- publicationGovernanceCopy(spec)
    changed$safety_limits$baseline_governors <- list("performance")
    expect_error(
        publicationValidateMeasureSpec(changed),
        class = "multischolar_publication_measure_spec_error"
    )

    missing <- publicationGovernanceCopy(spec)
    missing$safety_limits$frequency_telemetry_required <- FALSE
    expect_error(
        publicationValidateMeasureSpec(missing),
        class = "multischolar_publication_measure_spec_error"
    )

    unknown <- publicationGovernanceCopy(spec)
    unknown$safety_limits$safety_mode <- "unknown"
    expect_error(
        publicationValidateMeasureSpec(unknown),
        class = "multischolar_publication_measure_spec_error"
    )
})

test_that("cache labels fail closed without privileged cold evidence", {
    spec <- publicationMeasureSpecFixture()
    cache <- publicationReadJson(publicationCacheContractPath())

    fraud <- publicationGovernanceCopy(spec$cache_evidence)
    fraud$stratum <- "cold"
    expect_error(
        publicationValidateCacheEvidence(fraud, cache),
        class = "multischolar_publication_cache_error"
    )

    unverified <- publicationGovernanceCopy(spec$cache_evidence)
    unverified$stratum <- "verified_cold"
    unverified$primary_comparison_eligible <- FALSE
    expect_error(
        publicationValidateCacheEvidence(unverified, cache),
        class = "multischolar_publication_cache_error"
    )

    uncontrolled <- publicationGovernanceCopy(spec$cache_evidence)
    uncontrolled$stratum <- "cold_uncontrolled"
    uncontrolled$pre_read_complete <- FALSE
    uncontrolled$primary_comparison_eligible <- FALSE
    expect_silent(publicationValidateCacheEvidence(uncontrolled, cache))
})
