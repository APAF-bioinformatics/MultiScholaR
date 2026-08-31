publicationWriteCgroupFixture <- function(path, swap_bytes = 0) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    writeLines("104857600", file.path(path, "memory.current"))
    writeLines("125829120", file.path(path, "memory.peak"))
    writeLines(as.character(swap_bytes), file.path(path, "memory.swap.current"))
    writeLines(c(
        "anon 67108864",
        "file 25165824",
        "kernel 8388608",
        "slab 4194304"
    ), file.path(path, "memory.stat"))
    writeLines(c(
        "usage_usec 2500000",
        "user_usec 1500000",
        "system_usec 1000000"
    ), file.path(path, "cpu.stat"))
    writeLines(c(
        "low 0", "high 0", "max 0", "oom 0", "oom_kill 0",
        "oom_group_kill 0"
    ), file.path(path, "memory.events"))
    writeLines(c("populated 1", "frozen 0"), file.path(path, "cgroup.events"))
    writeLines(
        "8:0 rbytes=4096 wbytes=8192 rios=2 wios=3 dbytes=0 dios=0",
        file.path(path, "io.stat")
    )
    writeLines("some avg10=0.00 avg60=0.00 avg300=0.00 total=0", file.path(
        path,
        "memory.pressure"
    ))
    writeLines("", file.path(path, "cpu.pressure"))
    writeLines("", file.path(path, "io.pressure"))
    writeLines(character(), file.path(path, "cgroup.procs"))
    invisible(path)
}

publicationMeasurementSampleFixture <- function(
    elapsed_seconds,
    oom = 0,
    process_count = 1,
    swap_bytes = 0
) {
    list(
        elapsed_seconds = elapsed_seconds,
        memory = list(
            current_bytes = 1000,
            peak_bytes = 1200,
            swap_current_bytes = swap_bytes,
            stat = list(anon = 700, file = 200, kernel = 100, slab = 50)
        ),
        cpu = list(usage_usec = 1000, user_usec = 700, system_usec = 300),
        io = list(rbytes = 0, wbytes = 0, rios = 0, wios = 0,
            dbytes = 0, dios = 0),
        pressure = list(memory = list(), cpu = list(), io = list()),
        events = list(
            memory = list(low = 0, high = 0, max = 0, oom = oom,
                oom_kill = oom, oom_group_kill = 0),
            cgroup = list(populated = 1, frozen = 0)
        ),
        smaps = list(available = TRUE, pss_bytes = 900,
            uss_bytes = 800, rss_bytes = 1100),
        resources = list(file_descriptors = 4, thread_tasks = 1,
            kernel_locks = 0),
        process_count = process_count,
        pids = list(123L),
        disk = list(
            logical_bytes = 100,
            allocated_bytes = 4096,
            file_count = 1,
            category_logical_bytes = list(
                diagnostics = 0, staging_snapshot = 0, duckdb_spill = 0,
                committed = 0, final = 0, harness = 100
            )
        )
    )
}

publicationProcessFixture <- function(exit_status = 0L) {
    list(get_exit_status = \() as.integer(exit_status))
}

test_that("cgroup snapshot keeps charged anonymous file kernel and swap separate", {
    root <- withr::local_tempdir(pattern = "publication-cgroup-")
    publicationWriteCgroupFixture(root)

    snapshot <- publicationCgroupSnapshot(root)

    expect_identical(snapshot$memory$current_bytes, 104857600)
    expect_identical(snapshot$memory$peak_bytes, 125829120)
    expect_identical(snapshot$memory$stat$anon, 67108864)
    expect_identical(snapshot$memory$stat$file, 25165824)
    expect_identical(snapshot$memory$stat$kernel, 8388608)
    expect_identical(snapshot$memory$stat$slab, 4194304)
    expect_identical(snapshot$memory$swap_current_bytes, 0)
    expect_identical(snapshot$cpu$usage_usec, 2500000)
    expect_identical(snapshot$events$memory$oom_kill, 0)
    expect_identical(snapshot$events$cgroup$populated, 1)
    expect_identical(snapshot$io$rbytes, 4096)
    expect_identical(snapshot$io$wbytes, 8192)
    expect_false(snapshot$smaps$available)
    expect_null(publicationCgroupSnapshot(file.path(root, "missing")))
})

test_that("smaps rollup reports current process PSS USS and RSS", {
    smaps <- publicationReadSmapsRollup(Sys.getpid())

    expect_true(smaps$available)
    expect_gt(smaps$rss_bytes, 0)
    expect_gt(smaps$pss_bytes, 0)
    expect_gt(smaps$uss_bytes, 0)
    expect_lte(smaps$pss_bytes, smaps$rss_bytes)
    expect_lte(smaps$uss_bytes, smaps$rss_bytes)

    missing <- publicationReadSmapsRollup(2147483647L)
    expect_false(missing$available)

    truncated <- publicationParseSmapsRollup(c(
        "Rss: 100 kB",
        "Pss: 90 kB",
        "Private_Clean: 80 kB"
    ), 123L)
    expect_false(truncated$available)
    expect_true(all(c("Private_Dirty", "Private_Hugetlb") %in%
        unlist(truncated$missing_fields)))
})

test_that("disk accounting separates logical allocated and category bytes", {
    root <- withr::local_tempdir(pattern = "publication-disk-")
    dir.create(file.path(root, "committed"), recursive = TRUE)
    dir.create(file.path(root, "duckdb-spill"), recursive = TRUE)
    writeBin(as.raw(rep(1L, 100L)), file.path(root, "committed", "data.bin"))
    writeBin(as.raw(rep(2L, 50L)), file.path(root, "duckdb-spill", "spill.bin"))

    metrics <- publicationDirectoryMetrics(root)

    expect_identical(metrics$logical_bytes, 150)
    expect_gte(metrics$allocated_bytes, metrics$logical_bytes)
    expect_identical(metrics$file_count, 2L)
    expect_identical(metrics$category_logical_bytes$committed, 100)
    expect_identical(metrics$category_logical_bytes$duckdb_spill, 50)
})

test_that("host envelope and fallback retain certification boundaries", {
    host <- publicationHostEnvelope()
    fallback <- publicationEngineeringFallback("owned_cgroup_unavailable")

    expect_identical(host$schema, "multischolar.omics_publication_host")
    expect_identical(host$cgroup$version, 2L)
    expect_true(all(c("cpu", "io", "memory", "pids") %in%
        unlist(host$cgroup$controllers)))
    expect_gt(host$memory$total_bytes, 0)
    expect_true(publicationScalarNumber(host$cpu$sockets, positive = TRUE))
    expect_true(publicationScalarNumber(host$cpu$physical_cores, positive = TRUE))
    expect_true(publicationScalarNumber(host$cpu$logical_cores, positive = TRUE))
    expect_true(is.list(host$cpu$affinity))
    expect_true(is.list(host$numa))
    expect_true(is.list(host$compiler))
    expect_true(is.list(host$blas))
    expect_false(fallback$publication_certifiable)
    expect_false(fallback$promotion_authority)
})

test_that("metric contract binds units provenance and non-additive memory", {
    contract <- publicationReadJson(publicationMetricContractPath())
    binding <- publicationMetricContractBinding()

    expect_silent(publicationValidateMetricContract(contract))
    expect_match(binding$sha256, "^[0-9a-f]{64}$")
    expect_false(contract$rules$absent_is_zero)
    expect_true(contract$rules$memory_stat_slab_is_subset_of_kernel)
    expect_identical(
        contract$rules$memory_accounted_sum,
        "anon + file + kernel"
    )

    drift <- publicationGovernanceCopy(contract)
    drift$rules$memory_stat_slab_is_subset_of_kernel <- FALSE
    expect_error(
        publicationValidateMetricContract(drift),
        class = "multischolar_publication_schema_error"
    )
})

test_that("host preflight fails closed when owned systemd or safety checks fail", {
    host <- publicationHostEnvelope()
    contract <- publicationReadJson(
        "tests/testdata/omics-performance/host-preflight-contract-v2.json"
    )
    expect_silent(publicationValidateHostPreflightContract(contract))
    host$cgroup$systemd_user_available <- TRUE
    host$memory$available_bytes <- 64 * 1024^3
    host$filesystem <- list(
        available = TRUE,
        writable = TRUE,
        available_bytes = 512 * 1024^3
    )
    host$load_average <- list(0, 0, 0)

    passed <- publicationEvaluateHostPreflight(
        host,
        contract,
        expected_peak_bytes = 8 * 1024^3,
        maximum_temporary_bytes = 16 * 1024^3
    )
    expect_true(passed$certified)

    swapped <- publicationEvaluateHostPreflight(
        host,
        contract,
        expected_peak_bytes = 8 * 1024^3,
        maximum_temporary_bytes = 16 * 1024^3,
        swap_activity_bytes = 1
    )
    expect_false(swapped$certified)
    expect_identical(swapped$failure_outcome, "host_not_publication_certified")

    unavailable <- host
    unavailable$cgroup$systemd_user_available <- FALSE
    expect_false(publicationEvaluateHostPreflight(
        unavailable,
        contract,
        expected_peak_bytes = 8 * 1024^3,
        maximum_temporary_bytes = 16 * 1024^3
    )$certified)

    unknown_load <- host
    unknown_load$load_average <- list(NA_real_, 0, 0)
    expect_false(publicationEvaluateHostPreflight(
        unknown_load,
        contract,
        expected_peak_bytes = 8 * 1024^3,
        maximum_temporary_bytes = 16 * 1024^3
    )$certified)

    unknown_thermal <- host
    unknown_thermal$thermal <- list(available = FALSE)
    expect_false(publicationEvaluateHostPreflight(
        unknown_thermal,
        contract,
        expected_peak_bytes = 8 * 1024^3,
        maximum_temporary_bytes = 16 * 1024^3
    )$certified)

    unknown_frequency <- host
    unknown_frequency$cpu$frequencies <- list(available = FALSE)
    expect_false(publicationEvaluateHostPreflight(
        unknown_frequency,
        contract,
        expected_peak_bytes = 8 * 1024^3,
        maximum_temporary_bytes = 16 * 1024^3
    )$certified)

    contract_drift <- publicationGovernanceCopy(contract)
    contract_drift$required_controllers <- list("memory")
    expect_error(
        publicationValidateHostPreflightContract(contract_drift),
        class = "multischolar_publication_host_preflight_error"
    )
})

test_that("runtime safety detects swap disk memory load thermal and frequency", {
    sample <- list(
        memory = list(swap_current_bytes = 0),
        disk = list(allocated_bytes = 1024)
    )
    host <- list(
        memory = list(MemAvailable = 32 * 1024^3),
        filesystem = list(available = TRUE, available_bytes = 512 * 1024^3),
        load_average = list(0.1, 0.1, 0.1),
        frequencies = list(median_hz = 3000000000),
        thermal = list(available = TRUE, maximum_celsius = 60)
    )
    limits <- list(
        minimum_available_memory_bytes = 8 * 1024^3,
        minimum_available_disk_bytes = 100 * 1024^3,
        maximum_load = 1.6,
        maximum_thermal_celsius = 85,
        baseline_cpu_frequency_hz = 3000000000,
        maximum_cpu_frequency_drift_fraction = 0.05,
        maximum_run_allocated_disk_bytes = 50 * 1024^3
    )

    expect_true(publicationEvaluateRuntimeSafety(sample, host, limits)$safe)

    swapped <- sample
    swapped$memory$swap_current_bytes <- 1
    expect_false(publicationEvaluateRuntimeSafety(swapped, host, limits)$safe)

    loaded <- host
    loaded$load_average <- list(2, 2, 2)
    expect_false(publicationEvaluateRuntimeSafety(sample, loaded, limits)$safe)

    hot <- host
    hot$thermal$maximum_celsius <- 90
    expect_false(publicationEvaluateRuntimeSafety(sample, hot, limits)$safe)

    drifted <- host
    drifted$frequencies$median_hz <- 2000000000
    expect_false(publicationEvaluateRuntimeSafety(sample, drifted, limits)$safe)
})

test_that("dynamic host successor observes frequency and freezes governors", {
    expect_identical(
        publicationFileDigest(
            "tests/testdata/omics-performance/host-preflight-contract-v2.json"
        ),
        "4d95546438091bc3d0333443edb60135823c2b8d21bbb68834d6ab3273601c19"
    )
    contract <- publicationReadJson(
        "tests/testdata/omics-performance/host-preflight-contract-v3.json"
    )
    expect_silent(publicationValidateHostPreflightContract(contract))
    host <- publicationHostEnvelope()
    host$cgroup$systemd_user_available <- TRUE
    host$memory$available_bytes <- 64 * 1024^3
    host$filesystem <- list(
        available = TRUE,
        writable = TRUE,
        available_bytes = 512 * 1024^3
    )
    host$load_average <- list(0, 0, 0)
    host$cpu$frequencies <- list(
        available = TRUE,
        minimum_hz = 600000000,
        median_hz = 1200000000,
        maximum_hz = 3600000000
    )
    host$cpu$governors <- list("powersave")
    preflight <- publicationEvaluateHostPreflight(
        host,
        contract,
        expected_peak_bytes = 8 * 1024^3,
        maximum_temporary_bytes = 16 * 1024^3
    )
    expect_true(preflight$certified)
    expect_identical(preflight$power_policy$mode, "observe_do_not_modify")
    expect_false(preflight$power_policy$mutation_allowed)

    historical_load_host <- host
    historical_load_host$load_average <- list(0.1, 8, 8)
    expect_true(publicationEvaluateHostPreflight(
        historical_load_host,
        contract,
        expected_peak_bytes = 8 * 1024^3,
        maximum_temporary_bytes = 16 * 1024^3
    )$certified)

    limits <- publicationDynamicRuntimeSafetyLimits(
        host,
        preflight,
        maximum_run_allocated_disk_bytes = 50 * 1024^3
    )
    sample <- publicationMeasurementSampleFixture(1)
    loaded <- list(
        memory = list(MemAvailable = 32 * 1024^3),
        filesystem = list(available = TRUE, available_bytes = 512 * 1024^3),
        load_average = list(0.1, 0.1, 0.1),
        frequencies = list(
            available = TRUE,
            minimum_hz = 1800000000,
            median_hz = 3200000000,
            maximum_hz = 4100000000
        ),
        governors = list("powersave"),
        thermal = list(available = TRUE, maximum_celsius = 60)
    )
    observed <- publicationEvaluateRuntimeSafety(sample, loaded, limits)
    expect_true(observed$safe)
    expect_identical(observed$safety_mode, "dynamic_observed_v1")

    historical_load <- loaded
    historical_load$load_average <- list(0.1, 8, 8)
    expect_true(publicationEvaluateRuntimeSafety(
        sample,
        historical_load,
        limits
    )$safe)

    changed <- loaded
    changed$governors <- list("performance")
    expect_false(publicationEvaluateRuntimeSafety(sample, changed, limits)$safe)
    missing <- loaded
    missing$frequencies$available <- FALSE
    expect_false(publicationEvaluateRuntimeSafety(sample, missing, limits)$safe)
})

test_that("safety traces are bounded digest-bound and freshness-scoped", {
    record <- list(
        safe = TRUE,
        safety_mode = "dynamic_observed_v1",
        checks = list(frequency_telemetry = TRUE, governor_stability = TRUE),
        frequency_drift_fraction = NULL,
        observation = list(
            governors = list("powersave"),
            frequencies = list(
                available = TRUE,
                minimum_hz = 600000000,
                median_hz = 2400000000,
                maximum_hz = 4000000000
            ),
            load_average = list(0, 0, 0),
            maximum_thermal_celsius = 60,
            available_memory_bytes = 32 * 1024^3,
            available_disk_bytes = 512 * 1024^3
        ),
        reason = NULL,
        elapsed_seconds = 0
    )
    trace <- list(record, within(record, elapsed_seconds <- 1))
    first <- publicationSafetyTraceEvidence(trace, 10, 1)
    second <- publicationSafetyTraceEvidence(trace, 10, 1)
    expect_true(first$valid)
    expect_identical(first$record_count, 2L)
    expect_identical(first$trace_digest, second$trace_digest)
    expect_silent(publicationValidateSafetyTraceEvidence(first))
    changed <- trace
    changed[[2L]]$observation$frequencies$median_hz <- 3000000000
    expect_false(identical(
        first$trace_digest,
        publicationSafetyTraceEvidence(changed, 10, 1)$trace_digest
    ))
    tampered <- first
    tampered$trace[[2L]]$observation$frequencies$median_hz <- 3000000000
    expect_error(
        publicationValidateSafetyTraceEvidence(tampered),
        class = "multischolar_publication_host_preflight_error"
    )
    expect_false(publicationSafetyTraceEvidence(rep(trace, 6L), 10, 1)$valid)

    governors <- publicationCpuGovernors()
    skip_if_not(length(governors) > 0L)
    monitor <- publicationRuntimeSafetyMonitor(
        list(
            safety_mode = "dynamic_observed_v1",
            minimum_available_memory_bytes = 1,
            minimum_available_disk_bytes = 1,
            maximum_load = 1000,
            maximum_thermal_celsius = 1000,
            baseline_governors = governors,
            frequency_telemetry_required = TRUE,
            governor_stability_required = TRUE,
            maximum_run_allocated_disk_bytes = 50 * 1024^3
        ),
        interval_seconds = 1
    )
    sample <- publicationMeasurementSampleFixture(0)
    expect_true(monitor(sample, 0)$fresh)
    expect_false(monitor(sample, 0.5)$fresh)
    expect_true(monitor(sample, 1)$fresh)
})

test_that("pilot terminal states never decide from invalid measurements", {
    valid <- list(
        status = "passed",
        publication_certifiable = TRUE,
        timed_out = FALSE,
        safety_aborted = FALSE,
        phase_evidence = list(valid = TRUE),
        safety_evidence = publicationTestSafetyEvidence(),
        metrics = list(
            peak_charged_memory_bytes = 4 * 1024^3,
            elapsed_seconds = 30
        )
    )
    qualified <- publicationPilotQualificationDecision(valid)
    expect_true(qualified$decision_valid)
    expect_true(qualified$dimension_selection_allowed)
    expect_identical(
        publicationPilotTerminalStatus(valid, qualified),
        "qualified"
    )

    cases <- list(
        safety_aborted_no_dimension_decision = within(valid, {
            status <- "failed"
            publication_certifiable <- FALSE
            safety_aborted <- TRUE
        }),
        timed_out_no_dimension_decision = within(valid, {
            status <- "failed"
            publication_certifiable <- FALSE
            timed_out <- TRUE
        }),
        phase_invalid_no_dimension_decision = within(valid, {
            status <- "failed"
            publication_certifiable <- FALSE
            phase_evidence <- list(valid = FALSE)
        }),
        safety_evidence_invalid_no_dimension_decision = within(valid, {
            status <- "failed"
            publication_certifiable <- FALSE
            safety_evidence <- list(valid = FALSE)
        }),
        worker_failed_no_dimension_decision = within(valid, {
            status <- "failed"
            publication_certifiable <- FALSE
        })
    )
    for (expected in names(cases)) {
        decision <- publicationPilotQualificationDecision(cases[[expected]])
        expect_false(decision$decision_valid, info = expected)
        expect_false(decision$dimension_selection_allowed, info = expected)
        expect_identical(
            publicationPilotTerminalStatus(cases[[expected]], decision),
            expected
        )
    }
})

test_that("publication tooling cannot mutate host power policy", {
    files <- list.files(
        "tools/profiling",
        pattern = "[.]R$",
        full.names = TRUE
    )
    source <- unlist(lapply(files, readLines, warn = FALSE), use.names = FALSE)
    prohibited <- c("cpupower", "frequency-set", "sudo")
    expect_false(any(vapply(prohibited, \(value) {
        any(grepl(value, source, fixed = TRUE))
    }, logical(1))))
})

test_that("dynamic frequency successor binds every shared consumer", {
    record <- publicationReadJson(
        "tests/testdata/omics-performance/dynamic-frequency-successor-v2.json"
    )
    expect_identical(
        record$schema,
        "multischolar.omics_publication_dynamic_frequency_successor"
    )
    expect_identical(record$semantics$power_policy, "observe_do_not_modify")
    expect_false(record$semantics$power_policy_mutation_allowed)
    expect_false(record$semantics$idle_to_load_frequency_drift_aborts)
    expect_false(record$semantics$aborted_pilot_dimension_decision_allowed)
    expect_false(record$performance_authority)
    expect_false(record$promotion_authority)
    expect_false(record$publication_authority)
    expect_identical(
        publicationFileDigest(record$predecessor$path),
        record$predecessor$sha256
    )
    expect_identical(
        publicationFileDigest(record$successor$path),
        record$successor$sha256
    )
    for (binding in record$implementation_bindings) {
        expect_identical(
            publicationFileDigest(binding$path),
            binding$sha256,
            info = binding$path
        )
    }
    for (transition in record$kernel_transitions) {
        expect_true(publicationBindingMatchesCurrentOrSuccessor(list(
            path = transition$path,
            sha256 = transition$predecessor_sha256
        )), info = transition$path)
    }
    expect_true(all(vapply(record$excluded_attempts, \(attempt) {
        !isTRUE(attempt$dimension_decision_allowed)
    }, logical(1))))
})

test_that("retained windows exclude post-dwell work and retain trace provenance", {
    samples <- lapply(seq(2.9, 5.1, by = 0.05), publicationMeasurementSampleFixture)
    retained <- publicationRetainedSamples(samples, 0, 5, 2)
    diagnostics <- publicationRetainedDiagnostics(retained)
    window <- publicationRetainedWindowEvidence(
        retained,
        samples,
        0,
        5,
        2,
        0.1,
        0.5
    )

    expect_true(window$valid)
    expect_gte(min(vapply(retained, `[[`, numeric(1), "elapsed_seconds")), 3)
    expect_lte(max(vapply(retained, `[[`, numeric(1), "elapsed_seconds")), 5)
    expect_equal(diagnostics$charged_memory_slope_bytes_per_second, 0)
    expect_match(diagnostics$trace_sha256, "^[0-9a-f]{64}$")
    expect_length(diagnostics$trace, length(retained))

    jittered <- retained
    jittered[[length(jittered)]]$cpu$usage_usec <-
        jittered[[1L]]$cpu$usage_usec + 6
    expect_false(
        publicationRetainedDiagnostics(jittered)$background_activity_observed
    )
    jittered[[length(jittered)]]$cpu$usage_usec <-
        jittered[[1L]]$cpu$usage_usec + 2000
    expect_true(
        publicationRetainedDiagnostics(jittered)$background_activity_observed
    )
})

test_that("incomplete cgroup teardown snapshots are not sampled", {
    path <- withr::local_tempdir()

    expect_null(publicationSampleCgroup(
        path,
        elapsed_seconds = 1,
        disk = list(logical_bytes = 0, allocated_bytes = 0, file_count = 0)
    ))
})

test_that("PID ancestry distinguishes owned workers from injected processes", {
    own <- publicationPidOwnership(Sys.getpid(), Sys.getpid())
    injected <- publicationPidOwnership(Sys.getpid(), 2147483647L)
    empty <- publicationPidOwnership(integer(), Sys.getpid())

    expect_true(own$valid)
    expect_length(own$unexpected_pids, 0L)
    expect_length(own$indeterminate_pids, 0L)
    expect_false(injected$valid)
    expect_identical(injected$unexpected_pids, list(Sys.getpid()))
    expect_true(empty$valid)
    expect_length(empty$unexpected_pids, 0L)
})

test_that("timeout crash cgroup loss and OOM cannot remain certifiable", {
    samples <- lapply(seq(0, 5, by = 0.05), publicationMeasurementSampleFixture)
    execution <- list(
        sampling_interval_ms = 20,
        retained_dwell_seconds = 5,
        retained_window_seconds = 2,
        maximum_boundary_bracket_seconds = 0.5,
        retained_boundary_tolerance_ms = 100
    )
    state <- list(
        active_tasks = 0, open_queries = 0, open_writers = 0,
        open_leases = 0, open_connections = 0, open_results = 0,
        active_child_processes = 0,
        arrow_pool_bytes = 0,
        duckdb_memory_bytes = 0,
        duckdb_spill_bytes = 0,
        duckdb_prepared_statements = 0,
        temporary_paths = 0,
        cache_entries = 0,
        observers = 0,
        native_resources = 0,
        retained_dwell_seconds = 5,
        retention_acknowledgement = "fifo_v1",
        settled_monotonic_seconds = 100
    )
    empty_resources <- publicationEmptyWorkerResources()
    resource_ledger <- list(
        schema = "multischolar.omics_publication_worker_resources",
        schema_version = "1.0.0",
        high_water = empty_resources,
        retained = empty_resources,
        terminal = empty_resources
    )
    final_disk <- list(logical_bytes = 100, allocated_bytes = 4096, file_count = 1)
    measure <- \(process = publicationProcessFixture(), timed_out = FALSE,
        cgroup_lost = FALSE, values = samples) publicationCgroupResult(
        process = process,
        samples = values,
        marker_elapsed = 0,
        retention_state = state,
        execution = execution,
        timed_out = timed_out,
        cgroup_observed = TRUE,
        cgroup_lost = cgroup_lost,
        expected_root_pid = 123L,
        pid_ownership_valid = TRUE,
        unexpected_pids = integer(),
        indeterminate_pids = integer(),
        retention_acknowledged = TRUE,
        unit_result = list(
            available = TRUE,
            source = "systemd_unit_properties",
            result = "success",
            exec_main_code = 1L,
            exec_main_status = 0L,
            memory_peak_bytes = 1200,
            oom_killed = FALSE
        ),
        worker_resource_ledger = resource_ledger,
        cleanup = list(
            valid = TRUE,
            unit_active = FALSE,
            cgroup_exists = FALSE,
            elapsed_seconds = 0.01
        ),
        unit_name = "fixture",
        cgroup_path = "/fixture",
        stdout_path = "stdout",
        stderr_path = "stderr",
        started_at = "2026-08-27T00:00:00Z",
        total_elapsed = 5,
        final_disk = final_disk,
        safety_aborted = FALSE,
        safety_reason = NULL
    )

    passed <- measure()
    crashed <- measure(process = publicationProcessFixture(1L))
    timeout <- measure(timed_out = TRUE)
    lost <- measure(cgroup_lost = TRUE)
    oom_samples <- lapply(seq(0, 5, by = 0.05),
        \(elapsed) publicationMeasurementSampleFixture(elapsed, oom = 1))
    oom <- measure(values = oom_samples)
    child_samples <- lapply(seq(0, 5, by = 0.05),
        \(elapsed) publicationMeasurementSampleFixture(
            elapsed,
            process_count = 2
        ))
    child_leak <- measure(values = child_samples)
    nonfinite_samples <- publicationGovernanceCopy(samples)
    for (index in seq_along(nonfinite_samples)) {
        nonfinite_samples[[index]]$smaps$pss_bytes <- NA_real_
    }
    nonfinite <- measure(values = nonfinite_samples)
        background_samples <- publicationGovernanceCopy(samples)
        for (index in seq_along(background_samples)) {
            background_samples[[index]]$cpu$usage_usec <- 1000 + index * 100L
        }
    background <- measure(values = background_samples)
    swapped_samples <- lapply(seq(0, 5, by = 0.05),
        \(elapsed) publicationMeasurementSampleFixture(
            elapsed,
            swap_bytes = 4096
        ))
    swapped <- measure(values = swapped_samples)

    expect_identical(passed$status, "passed")
    expect_true(passed$publication_certifiable)
    for (result in list(
        crashed, timeout, lost, oom, child_leak, nonfinite, swapped, background
    )) {
        expect_identical(result$status, "failed")
        expect_false(result$publication_certifiable)
    }

    leaked_resources <- publicationGovernanceCopy(resource_ledger)
    leaked_resources$terminal$duckdb_connections <- 1
    expect_false(publicationWorkerResourceLedgerValid(leaked_resources))
})

test_that("owned user cgroup integration measures null worker exactly", {
    skip_if_not(identical(
        Sys.getenv("MULTISCHOLAR_PUBLICATION_CGROUP_TEST"),
        "true"
    ))
    skip_if_not(publicationSystemdUserAvailable())

    output <- tempfile("null-calibration-", fileext = ".json")
    process <- processx::run(
        file.path(R.home("bin"), "Rscript"),
        c(
            "--vanilla",
            publicationToolPath("run_omics_publication_benchmark.R"),
            "--null-calibration", "true",
            "--runs", "2",
            "--dwell-seconds", "5",
            "--output", output
        ),
        wd = publicationGovernanceRepoPath(),
        timeout = 120000,
        error_on_status = FALSE,
        echo = FALSE
    )
    expect_identical(
        process$status,
        0L,
        info = paste("stdout:", process$stdout, "stderr:", process$stderr)
    )

    result <- jsonlite::read_json(output, simplifyVector = FALSE)
    failure_evidence <- paste(vapply(result$runs, \(run) {
        jsonlite::toJSON(list(
            run_id = run$run_id,
            status = run$status,
            exit_status = run$exit_status,
            timed_out = run$timed_out,
            safety_aborted = run$safety_aborted,
            safety_reason = run$safety_reason,
            cgroup_observed = run$cgroup_observed,
            cgroup_lost = run$cgroup_lost,
            pid_ownership_valid = run$pid_ownership_valid,
            unexpected_pids = run$unexpected_pids,
            indeterminate_pids = run$indeterminate_pids,
            unit_result = run$unit_result,
            worker_resource_ledger_valid = run$worker_resource_ledger_valid,
            cleanup = run$cleanup,
            retention_state_valid = run$retention_state_valid,
            retained_window = run$retained_window,
            retained_process_count = run$retained_diagnostics$maximum_process_count,
            metrics_valid = run$metrics_valid,
            memory_events = run$metrics$memory_events
        ), auto_unbox = TRUE, null = "null")
    }, character(1)), collapse = "\n")
    expect_identical(result$status, "passed", info = failure_evidence)
    expect_length(result$runs, 2L)
    expect_true(all(vapply(result$runs, \(run) {
        identical(run$status, "passed") &&
            isTRUE(run$cgroup_observed) &&
            isTRUE(run$publication_certifiable) &&
            isTRUE(run$retention_state_valid) &&
            isTRUE(run$worker_resource_ledger_valid) &&
            isTRUE(run$metrics_valid) &&
            isTRUE(run$phase_evidence$valid) &&
            run$metrics$peak_charged_memory_bytes > 0 &&
            run$metrics$retained_charged_memory_bytes > 0 &&
            run$metrics$peak_anonymous_memory_bytes > 0 &&
            run$metrics$peak_pss_bytes > 0 &&
            run$metrics$peak_uss_bytes > 0 &&
            run$metrics$peak_rss_bytes > 0 &&
            run$metrics$phase_cpu_seconds > 0 &&
            run$metrics$cpu_usage_seconds > 0 &&
            run$metrics$primary_work_units_per_wall_second > 0 &&
            run$metrics$primary_work_units_per_cpu_second > 0 &&
            run$metrics$peak_file_descriptors > 0 &&
            run$metrics$retained_file_descriptors > 0 &&
            run$metrics$peak_thread_tasks > 0 &&
            run$metrics$retained_thread_tasks > 0 &&
            run$metrics$peak_kernel_locks >= 0 &&
            run$metrics$peak_logical_disk_bytes > 0 &&
            run$metrics$peak_allocated_disk_bytes > 0 &&
            run$metrics$final_logical_disk_bytes > 0 &&
            run$metrics$final_allocated_disk_bytes > 0 &&
            run$metrics$final_file_count > 0 &&
            run$metrics$peak_swap_bytes == 0
    }, logical(1))), info = failure_evidence)
})

test_that("owned cgroup terminates unsafe workload without certification", {
    skip_if_not(identical(
        Sys.getenv("MULTISCHOLAR_PUBLICATION_CGROUP_TEST"),
        "true"
    ))
    skip_if_not(publicationSystemdUserAvailable())

    run_dir <- tempfile("unsafe-cgroup-")
    result <- publicationMeasureCgroupSubprocess(
        command = file.path(R.home("bin"), "Rscript"),
        command_args = c(
            "--vanilla",
            publicationToolPath("run_omics_publication_benchmark.R"),
            "--worker", "null",
            "--run-dir", run_dir,
            "--dwell-seconds", "5"
        ),
        run_dir = run_dir,
        execution = list(
            sampling_interval_ms = 20,
            disk_sampling_interval_ms = 500,
            timeout_seconds = 30,
            retained_dwell_seconds = 5,
            retained_window_seconds = 2
        ),
        safety_check_fn = \(...) list(safe = FALSE, reason = "injected_unsafe"),
        unit_name = publicationSystemdUnitName("multischolar-unsafe")
    )

    expect_identical(result$status, "failed")
    expect_true(result$safety_aborted)
    expect_identical(result$safety_reason, "injected_unsafe")
    expect_false(result$publication_certifiable)
})

test_that("owned cgroup records crash timeout child leak and bounded OOM", {
    skip_if_not(identical(
        Sys.getenv("MULTISCHOLAR_PUBLICATION_CGROUP_TEST"),
        "true"
    ))
    skip_if_not(publicationSystemdUserAvailable())

    run_fault <- \(mode, timeout_seconds = 20, memory_max_bytes = NULL) {
        run_dir <- tempfile(paste0("fault-", mode, "-"))
        execution <- list(
            sampling_interval_ms = 20,
            disk_sampling_interval_ms = 500,
            timeout_seconds = timeout_seconds,
            retained_dwell_seconds = 5,
            retained_window_seconds = 2,
            maximum_boundary_bracket_seconds = 0.5,
            retained_boundary_tolerance_ms = 100
        )
        if (!is.null(memory_max_bytes)) {
            execution$memory_max_bytes <- memory_max_bytes
        }
        publicationMeasureCgroupSubprocess(
            command = file.path(R.home("bin"), "Rscript"),
            command_args = c(
                "--vanilla",
                publicationToolPath("run_omics_publication_benchmark.R"),
                "--worker", mode,
                "--run-dir", run_dir,
                "--dwell-seconds", "5"
            ),
            run_dir = run_dir,
            execution = execution,
            env = c(MULTISCHOLAR_PUBLICATION_FAULT_INJECTION = "true"),
            unit_name = publicationSystemdUnitName(paste0(
                "multischolar-fault-",
                mode
            ))
        )
    }

    crashed <- run_fault("crash")
    timed_out <- run_fault("timeout", timeout_seconds = 1)
    child_leak <- run_fault("child_leak")
    oom_binary <- tempfile("omics-publication-oom-probe-")
    compile <- processx::run(
        Sys.getenv("CC", unset = "cc"),
        c(
            "-O0", "-std=c11",
            publicationToolPath("omics_publication_oom_probe.c"),
            "-o", oom_binary
        ),
        error_on_status = FALSE,
        echo = FALSE
    )
    expect_identical(
        compile$status,
        0L,
        info = paste("stdout:", compile$stdout, "stderr:", compile$stderr)
    )
    oom_dir <- tempfile("fault-oom-")
    oom <- publicationMeasureCgroupSubprocess(
        command = oom_binary,
        command_args = character(),
        run_dir = oom_dir,
        execution = list(
            sampling_interval_ms = 20,
            disk_sampling_interval_ms = 500,
            timeout_seconds = 20,
            retained_dwell_seconds = 5,
            retained_window_seconds = 2,
            maximum_boundary_bracket_seconds = 0.5,
            retained_boundary_tolerance_ms = 100,
            memory_max_bytes = 128 * 1024^2
        ),
        unit_name = publicationSystemdUnitName("multischolar-fault-oom")
    )

    for (result in list(crashed, timed_out, child_leak, oom)) {
        expect_identical(result$status, "failed")
        expect_false(result$publication_certifiable)
        expect_true(result$cleanup$valid)
    }
    expect_false(identical(crashed$exit_status, 0L))
    expect_true(timed_out$timed_out)
    expect_gte(child_leak$retained_diagnostics$maximum_process_count, 2)
    expect_true(oom$unit_result$oom_killed)
    expect_identical(oom$unit_result$result, "oom-kill")
    expect_identical(oom$unit_result$exec_main_code, 2L)
    expect_identical(oom$unit_result$exec_main_status, 9L)
    expect_equal(oom$unit_result$memory_peak_bytes, 128 * 1024^2)
})

test_that("owned cgroup separates retained mmap page cache from anonymous memory", {
    skip_if_not(identical(
        Sys.getenv("MULTISCHOLAR_PUBLICATION_CGROUP_TEST"),
        "true"
    ))
    skip_if_not(publicationSystemdUserAvailable())

    binary <- tempfile("omics-publication-mmap-probe-")
    compile <- processx::run(
        Sys.getenv("CC", unset = "cc"),
        c(
            "-O2", "-std=c11",
            testthat::test_path("fixtures", "omics_publication_mmap_probe.c"),
            "-o", binary
        ),
        error_on_status = FALSE,
        echo = FALSE
    )
    expect_identical(
        compile$status,
        0L,
        info = paste("stdout:", compile$stdout, "stderr:", compile$stderr)
    )
    run_dir <- tempfile("mmap-page-cache-")
    result <- publicationMeasureCgroupSubprocess(
        command = binary,
        command_args = run_dir,
        run_dir = run_dir,
        execution = list(
            sampling_interval_ms = 20,
            disk_sampling_interval_ms = 500,
            timeout_seconds = 30,
            retained_dwell_seconds = 5,
            retained_window_seconds = 2,
            maximum_boundary_bracket_seconds = 0.5,
            retained_boundary_tolerance_ms = 100
        ),
        unit_name = publicationSystemdUnitName("multischolar-mmap-probe")
    )

    expect_identical(result$status, "passed")
    expect_true(result$publication_certifiable)
    expect_true(result$retention_acknowledged)
    expect_gt(result$metrics$retained_file_memory_bytes, 16 * 1024^2)
    expect_gt(result$metrics$retained_anonymous_memory_bytes, 0)
    expect_gt(result$metrics$retained_kernel_memory_bytes, 0)
    expect_gte(
        result$metrics$retained_kernel_memory_bytes,
        result$metrics$retained_slab_memory_bytes
    )
    expect_gt(result$metrics$retained_pss_bytes, 0)
    expect_gt(result$metrics$retained_uss_bytes, 0)
    expect_gt(result$metrics$retained_rss_bytes, 0)
    expect_identical(result$metrics$peak_swap_bytes, 0)
    expect_false(result$retained_diagnostics$background_activity_observed)
})
