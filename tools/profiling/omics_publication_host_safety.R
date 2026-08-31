publicationCpuModel <- function() {
    line <- grep(
        "^model name[[:space:]]*:",
        publicationReadLinesSafe("/proc/cpuinfo"),
        value = TRUE
    )
    if (!length(line)) return(NULL)
    trimws(sub("^[^:]+:", "", line[[1L]]))
}

publicationCpuTopology <- function() {
    paths <- Sys.glob("/sys/devices/system/cpu/cpu[0-9]*/topology/core_id")
    records <- lapply(paths, \(path) {
        cpu_path <- dirname(dirname(path))
        list(
            cpu = basename(cpu_path),
            core_id = publicationReadScalarFile(path),
            package_id = publicationReadScalarFile(file.path(
                dirname(path),
                "physical_package_id"
            ))
        )
    })
    valid <- vapply(records, \(record) {
        is.finite(record$core_id) && is.finite(record$package_id)
    }, logical(1))
    records <- records[valid]
    core_keys <- unique(vapply(records, \(record) {
        paste(record$package_id, record$core_id, sep = ":")
    }, character(1)))
    packages <- unique(vapply(records, `[[`, numeric(1), "package_id"))
    microcode <- grep(
        "^microcode[[:space:]]*:",
        publicationReadLinesSafe("/proc/cpuinfo"),
        value = TRUE
    )
    list(
        topology_available = length(records) > 0L,
        sockets = if (length(packages)) length(packages) else NA_integer_,
        physical_cores = if (length(core_keys)) length(core_keys) else NA_integer_,
        logical_cores = length(paths),
        microcode = as.list(unique(trimws(sub("^[^:]+:", "", microcode))))
    )
}

publicationCpuAffinity <- function() {
    lines <- publicationReadLinesSafe("/proc/self/status")
    value <- \(name) {
        match <- grep(paste0("^", name, ":[[:space:]]*"), lines, value = TRUE)
        if (length(match) != 1L) return(NULL)
        trimws(sub("^[^:]+:", "", match[[1L]]))
    }
    list(
        cpus_allowed_list = value("Cpus_allowed_list"),
        memory_nodes_allowed_list = value("Mems_allowed_list")
    )
}

publicationNumaTopology <- function() {
    paths <- Sys.glob("/sys/devices/system/node/node[0-9]*")
    nodes <- lapply(paths, \(path) {
        meminfo <- publicationReadLinesSafe(file.path(path, "meminfo"))
        total <- grep("MemTotal:", meminfo, value = TRUE)
        list(
            node = basename(path),
            cpu_list = paste(publicationReadLinesSafe(file.path(
                path,
                "cpulist"
            )), collapse = ""),
            memory_total_bytes = if (length(total)) {
                as.numeric(gsub("[^0-9]", "", total[[1L]])) * 1024
            } else {
                NA_real_
            }
        )
    })
    list(available = length(nodes) > 0L, node_count = length(nodes), nodes = nodes)
}

publicationRConfig <- function(name) {
    result <- processx::run(
        file.path(R.home("bin"), "R"),
        c("CMD", "config", name),
        error_on_status = FALSE,
        echo = FALSE
    )
    if (result$status != 0L) return(NULL)
    trimws(result$stdout)
}

publicationCompilerMetadata <- function() {
    compiler <- publicationRConfig("CC")
    command <- if (publicationScalarString(compiler)) {
        strsplit(compiler, "[[:space:]]+", perl = TRUE)[[1L]][[1L]]
    } else {
        NULL
    }
    version <- if (!is.null(command) && nzchar(Sys.which(command))) {
        result <- processx::run(
            command,
            "--version",
            error_on_status = FALSE,
            echo = FALSE
        )
        if (result$status == 0L) strsplit(result$stdout, "\n", fixed = TRUE)[[1L]][[1L]] else
            NULL
    } else {
        NULL
    }
    list(
        cc = compiler,
        cc_version = version,
        cxx = publicationRConfig("CXX"),
        cflags = publicationRConfig("CFLAGS"),
        cxxflags = publicationRConfig("CXXFLAGS")
    )
}

publicationBlasMetadata <- function() {
    software <- as.list(extSoftVersion())
    info <- sessionInfo()
    list(
        blas = info$BLAS,
        lapack = info$LAPACK,
        extsoft_blas = software$BLAS,
        extsoft_lapack = software$LAPACK
    )
}

publicationMeminfo <- function() {
    lines <- publicationReadLinesSafe("/proc/meminfo")
    values <- list()
    for (line in lines) {
        pair <- strsplit(line, ":", fixed = TRUE)[[1L]]
        if (length(pair) != 2L) next
        value <- suppressWarnings(as.numeric(gsub("[^0-9]", "", pair[[2L]])))
        if (is.finite(value)) values[[pair[[1L]]]] <- value * 1024
    }
    values
}

publicationCpuFrequencies <- function() {
    paths <- Sys.glob("/sys/devices/system/cpu/cpu*/cpufreq/scaling_cur_freq")
    values <- vapply(paths, publicationReadScalarFile, numeric(1)) * 1000
    values <- values[is.finite(values)]
    list(
        available = length(values) > 0L,
        minimum_hz = if (length(values)) min(values) else NA_real_,
        median_hz = if (length(values)) stats::median(values) else NA_real_,
        maximum_hz = if (length(values)) max(values) else NA_real_
    )
}

publicationCpuGovernors <- function() {
    paths <- Sys.glob("/sys/devices/system/cpu/cpu*/cpufreq/scaling_governor")
    values <- unique(unlist(lapply(paths, publicationReadLinesSafe), use.names = FALSE))
    as.list(values[nzchar(values)])
}

publicationNormalizeGovernors <- function(governors) {
    values <- as.character(unlist(governors, use.names = FALSE))
    sort(unique(values[!is.na(values) & nzchar(values)]), method = "radix")
}

publicationThermalState <- function() {
    paths <- Sys.glob("/sys/class/thermal/thermal_zone*/temp")
    values <- vapply(paths, publicationReadScalarFile, numeric(1)) / 1000
    values <- values[is.finite(values)]
    list(
        available = length(values) > 0L,
        minimum_celsius = if (length(values)) min(values) else NA_real_,
        median_celsius = if (length(values)) stats::median(values) else NA_real_,
        maximum_celsius = if (length(values)) max(values) else NA_real_
    )
}

publicationLoadAverage <- function() {
    lines <- publicationReadLinesSafe("/proc/loadavg")
    if (!length(lines)) return(as.list(rep(NA_real_, 3L)))
    fields <- strsplit(trimws(lines[[1L]]), "[[:space:]]+", perl = TRUE)[[1L]]
    values <- suppressWarnings(as.numeric(fields[seq_len(min(3L, length(fields)))]))
    if (length(values) < 3L) values <- c(values, rep(NA_real_, 3L - length(values)))
    as.list(values)
}

publicationFilesystemMetrics <- function(path = .PUBLICATION_REPO_ROOT) {
    output <- processx::run(
        "df",
        c("-P", "-B1", normalizePath(path, mustWork = TRUE)),
        error_on_status = FALSE,
        echo = FALSE
    )
    lines <- strsplit(trimws(output$stdout), "\n", fixed = TRUE)[[1L]]
    if (output$status != 0L || length(lines) < 2L) {
        return(list(available = FALSE))
    }
    fields <- strsplit(trimws(lines[[length(lines)]]), "[[:space:]]+", perl = TRUE)[[1L]]
    if (length(fields) < 6L) return(list(available = FALSE))
    list(
        available = TRUE,
        filesystem = fields[[1L]],
        total_bytes = as.numeric(fields[[2L]]),
        used_bytes = as.numeric(fields[[3L]]),
        available_bytes = as.numeric(fields[[4L]]),
        mount_point = fields[[6L]],
        writable = publicationFilesystemWritable(path)
    )
}

publicationFilesystemWritable <- function(path) {
    probe <- tempfile(".multischolar-write-probe-", tmpdir = path)
    created <- suppressWarnings(file.create(probe))
    if (isTRUE(created)) unlink(probe, force = TRUE)
    isTRUE(created)
}

publicationStorageMetadata <- function(path = .PUBLICATION_REPO_ROOT) {
    if (!nzchar(Sys.which("findmnt"))) return(list(available = FALSE))
    result <- processx::run(
        "findmnt",
        c(
            "--json", "--target", normalizePath(path, mustWork = TRUE),
            "--output", "SOURCE,FSTYPE,OPTIONS,SIZE,AVAIL,FSROOT,TARGET"
        ),
        error_on_status = FALSE,
        echo = FALSE
    )
    if (result$status != 0L) return(list(available = FALSE))
    parsed <- tryCatch(
        jsonlite::fromJSON(result$stdout, simplifyVector = FALSE),
        error = \(...) NULL
    )
    if (!is.list(parsed) || !length(parsed$filesystems)) {
        return(list(available = FALSE))
    }
    list(available = TRUE, filesystem = parsed$filesystems[[1L]])
}

publicationGitDirty <- function() {
    result <- processx::run(
        "git",
        c("status", "--porcelain", "--untracked-files=no"),
        wd = .PUBLICATION_REPO_ROOT,
        error_on_status = FALSE,
        echo = FALSE
    )
    if (result$status != 0L) return(NA)
    nzchar(trimws(result$stdout))
}

publicationSystemdUserAvailable <- function() {
    if (!nzchar(Sys.which("systemd-run")) || !nzchar(Sys.which("systemctl"))) {
        return(FALSE)
    }
    result <- processx::run(
        "systemctl",
        c("--user", "show-environment"),
        error_on_status = FALSE,
        echo = FALSE,
        timeout = 5000
    )
    identical(result$status, 0L)
}

publicationHostEnvelope <- function(filesystem_path = .PUBLICATION_REPO_ROOT) {
    info <- Sys.info()
    meminfo <- publicationMeminfo()
    topology <- publicationCpuTopology()
    cgroup_path <- publicationCurrentCgroupPath()
    controllers <- publicationReadLinesSafe(file.path(
        publicationCgroupRoot(),
        "cgroup.controllers"
    ))
    list(
        schema = "multischolar.omics_publication_host",
        schema_version = "1.0.0",
        captured_at = publicationUtcNow(),
        git_revision = publicationGitRevision(),
        git_dirty = publicationGitDirty(),
        r_version = R.version.string,
        platform = R.version$platform,
        os = as.list(info[intersect(
            c("sysname", "release", "version", "machine"),
            names(info)
        )]),
        cpu = list(
            model = publicationCpuModel(),
            sockets = topology$sockets,
            logical_cores = topology$logical_cores,
            physical_cores = topology$physical_cores,
            microcode = topology$microcode,
            affinity = publicationCpuAffinity(),
            frequencies = publicationCpuFrequencies(),
            governors = publicationCpuGovernors()
        ),
        numa = publicationNumaTopology(),
        memory = list(
            total_bytes = meminfo$MemTotal,
            available_bytes = meminfo$MemAvailable,
            swap_total_bytes = meminfo$SwapTotal,
            swap_free_bytes = meminfo$SwapFree
        ),
        cgroup = list(
            version = if (file.exists(file.path(
                publicationCgroupRoot(),
                "cgroup.controllers"
            ))) 2L else NA_integer_,
            current_path = cgroup_path,
            controllers = if (length(controllers)) {
                as.list(strsplit(controllers[[1L]], " ", fixed = TRUE)[[1L]])
            } else {
                list()
            },
            systemd_user_available = publicationSystemdUserAvailable()
        ),
        filesystem = publicationFilesystemMetrics(filesystem_path),
        storage = publicationStorageMetadata(filesystem_path),
        thermal = publicationThermalState(),
        load_average = publicationLoadAverage(),
        locale = Sys.getlocale(),
        timezone = Sys.timezone(),
        rng = stats::setNames(as.list(RNGkind()), c(
            "kind", "normal_kind", "sample_kind"
        )),
        compiler = publicationCompilerMetadata(),
        blas = publicationBlasMetadata(),
        thread_environment = as.list(Sys.getenv(
            c(
                "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                "ARROW_NUM_THREADS", "DUCKDB_THREADS"
            ),
            unset = NA_character_
        )),
        package_versions = publicationPackageVersions(c(
            "processx", "ps", "jsonlite", "digest", "arrow", "duckdb", "DBI"
        ))
    )
}

publicationValidateHostPreflightContract <- function(contract) {
    legacy_fields <- c(
        "schema", "schema_version", "preflight_contract_id", "owner_ticket_id",
        "platform", "required_cgroup_version", "required_controllers",
        "required_host_metadata", "owned_transient_unit", "preflight_rules",
        "mid_run_abort_rules", "certification", "unavailable_outcome"
    )
    successor <- identical(contract$schema_version, "1.1.0")
    expected <- if (successor) {
        append(legacy_fields, "predecessor", after = 4L)
    } else {
        legacy_fields
    }
    publicationRequireNames(contract, expected, "Host preflight contract")
    controllers <- unlist(contract$required_controllers, use.names = FALSE)
    metadata <- unlist(contract$required_host_metadata, use.names = FALSE)
    legacy_metadata <- c(
        "cpu_topology", "cpu_affinity", "numa_topology",
        "storage_filesystem", "compiler", "blas_lapack", "thermal",
        "cpu_frequency"
    )
    required_metadata <- if (successor) {
        c(legacy_metadata, "cpu_governor")
    } else {
        legacy_metadata
    }
    version_valid <- if (successor) {
        publicationRequireNames(contract$predecessor, c(
            "preflight_contract_id", "path", "sha256"
        ), "Host preflight predecessor")
        identical(contract$owner_ticket_id, "OMICS-ART-085") &&
            identical(
                contract$predecessor$preflight_contract_id,
                "multischolar.omics_publication_host_preflight.2026-08-27.v2"
            ) && identical(
                contract$predecessor$path,
                "tests/testdata/omics-performance/host-preflight-contract-v2.json"
            ) && identical(
                contract$predecessor$sha256,
                "4d95546438091bc3d0333443edb60135823c2b8d21bbb68834d6ab3273601c19"
            ) && identical(
                publicationFileDigest(contract$predecessor$path),
                contract$predecessor$sha256
            ) && identical(
                contract$preflight_rules$cpu_power_policy_rule,
                "observe_do_not_modify_and_require_stable_governor_set"
            ) && identical(
                contract$preflight_rules$cpu_frequency_rule,
                "telemetry_required_dynamic_frequency_is_observed_not_aborted"
            ) && identical(
                contract$preflight_rules$load_window_rule,
                "one_minute_current_window_with_longer_windows_recorded_only"
            ) && identical(
                contract$certification$power_policy_mutation_allowed,
                FALSE
            ) && isTRUE(
                contract$certification$governor_stability_required
            ) && isTRUE(
                contract$certification$frequency_telemetry_required
            )
    } else {
        identical(contract$schema_version, "1.0.0") &&
            identical(contract$owner_ticket_id, "OMICS-ART-063")
    }
    valid <- identical(
        contract$schema,
        "multischolar.omics_publication_host_preflight_contract"
    ) && isTRUE(version_valid) &&
        identical(contract$platform, "Linux") &&
        identical(contract$required_cgroup_version, 2L) &&
        setequal(controllers, c("cpu", "io", "memory", "pids")) &&
        setequal(metadata, required_metadata) &&
        isTRUE(contract$owned_transient_unit$memory_accounting) &&
        isTRUE(contract$owned_transient_unit$io_accounting) &&
        isTRUE(contract$owned_transient_unit$single_worker_per_unit) &&
        isTRUE(contract$certification$owned_cgroup_observed_required) &&
        isTRUE(contract$certification$five_second_retained_dwell_required) &&
        identical(
            contract$unavailable_outcome,
            "host_not_publication_certified"
        )
    if (!valid) {
        publicationAbort(
            "Host preflight contract differs",
            "multischolar_publication_host_preflight_error"
        )
    }
    invisible(contract)
}

publicationEvaluateHostPreflight <- function(
    host,
    contract,
    expected_peak_bytes,
    maximum_temporary_bytes,
    swap_activity_bytes = 0,
    cpu_frequency_drift_fraction = 0,
    thermal_throttled = FALSE,
    stale_owned_units = 0L,
    cleanup_receipt_current = TRUE
) {
    publicationValidateHostPreflightContract(contract)
    if (!publicationScalarNumber(expected_peak_bytes, positive = TRUE) ||
        !publicationScalarNumber(maximum_temporary_bytes, positive = TRUE)) {
        publicationAbort(
            "Host preflight size bounds are invalid",
            "multischolar_publication_host_preflight_error"
        )
    }
    required_controllers <- unlist(contract$required_controllers, use.names = FALSE)
    observed_controllers <- unlist(host$cgroup$controllers, use.names = FALSE)
    minimum_memory <- max(8589934592, 2 * expected_peak_bytes)
    minimum_disk <- max(107374182400, 2 * maximum_temporary_bytes)
    load_limit <- host$cpu$logical_cores *
        contract$preflight_rules$maximum_load_per_logical_core
    load_values <- as.numeric(unlist(host$load_average, use.names = FALSE))
    successor <- identical(contract$schema_version, "1.1.0")
    governors <- publicationNormalizeGovernors(host$cpu$governors)
    frequency_valid <- isTRUE(host$cpu$frequencies$available) &&
        publicationScalarNumber(
            host$cpu$frequencies$minimum_hz,
            positive = TRUE
        ) && publicationScalarNumber(
            host$cpu$frequencies$median_hz,
            positive = TRUE
        ) && publicationScalarNumber(
            host$cpu$frequencies$maximum_hz,
            positive = TRUE
        )
    frequency_check <- if (successor) {
        frequency_valid && length(governors) > 0L
    } else {
        frequency_valid && publicationScalarNumber(
            cpu_frequency_drift_fraction
        ) && cpu_frequency_drift_fraction <=
            contract$preflight_rules$maximum_cpu_frequency_drift_fraction
    }
    load_check <- length(load_values) == 3L && all(is.finite(load_values)) &&
        if (successor) {
            load_values[[1L]] <= load_limit
        } else {
            all(load_values <= load_limit)
        }
    checks <- list(
        linux = identical(host$os$sysname, "Linux"),
        cpu_topology = publicationScalarNumber(
            host$cpu$physical_cores,
            positive = TRUE
        ) && publicationScalarNumber(host$cpu$logical_cores, positive = TRUE) &&
            publicationScalarNumber(host$cpu$sockets, positive = TRUE),
        cpu_affinity = publicationScalarString(
            host$cpu$affinity$cpus_allowed_list
        ),
        numa = isTRUE(host$numa$available) &&
            publicationScalarNumber(host$numa$node_count, positive = TRUE),
        cgroup_v2 = identical(host$cgroup$version, 2L),
        controllers = all(required_controllers %in% observed_controllers),
        systemd_user = isTRUE(host$cgroup$systemd_user_available),
        memory = publicationScalarNumber(host$memory$available_bytes) &&
            host$memory$available_bytes >= minimum_memory,
        disk = isTRUE(host$filesystem$available) &&
            isTRUE(host$filesystem$writable) &&
            publicationScalarNumber(host$filesystem$available_bytes) &&
            host$filesystem$available_bytes >= minimum_disk,
        storage = isTRUE(host$storage$available),
        swap = publicationScalarNumber(swap_activity_bytes) &&
            swap_activity_bytes == 0,
        load = load_check,
        frequency = frequency_check,
        thermal = isTRUE(host$thermal$available) &&
            publicationScalarNumber(host$thermal$maximum_celsius) &&
            !isTRUE(thermal_throttled),
        compiler = publicationScalarString(host$compiler$cc) &&
            publicationScalarString(host$compiler$cc_version),
        blas = publicationScalarString(host$blas$blas) &&
            publicationScalarString(host$blas$lapack),
        stale_units = identical(as.integer(stale_owned_units), 0L),
        cleanup = isTRUE(cleanup_receipt_current)
    )
    result <- list(
        schema = "multischolar.omics_publication_host_preflight",
        schema_version = if (successor) "1.1.0" else "1.0.0",
        contract_id = contract$preflight_contract_id,
        checked_at = publicationUtcNow(),
        expected_peak_bytes = expected_peak_bytes,
        required_available_memory_bytes = minimum_memory,
        maximum_temporary_bytes = maximum_temporary_bytes,
        required_available_disk_bytes = minimum_disk,
        checks = checks,
        certified = isTRUE(all(unlist(checks, use.names = FALSE))),
        failure_outcome = if (isTRUE(all(unlist(checks, use.names = FALSE)))) NULL else
            contract$unavailable_outcome
    )
    if (successor) {
        result$power_policy <- list(
            mode = "observe_do_not_modify",
            governors = as.list(governors),
            frequency = host$cpu$frequencies,
            mutation_allowed = FALSE
        )
    }
    result
}

publicationEngineeringFallback <- function(reason, metrics = list()) {
    list(
        schema = "multischolar.omics_publication_engineering_fallback",
        schema_version = "1.0.0",
        evidence_class = "engineering_polling_only",
        reason = reason,
        metrics = metrics,
        publication_certifiable = FALSE,
        promotion_authority = FALSE
    )
}

publicationLightweightHostSnapshot <- function(
    filesystem_path = .PUBLICATION_REPO_ROOT
) {
    list(
        memory = publicationMeminfo(),
        filesystem = publicationFilesystemMetrics(filesystem_path),
        load_average = publicationLoadAverage(),
        frequencies = publicationCpuFrequencies(),
        governors = publicationCpuGovernors(),
        thermal = publicationThermalState()
    )
}

publicationDynamicRuntimeSafetyLimits <- function(
    host,
    preflight,
    maximum_run_allocated_disk_bytes
) {
    governors <- publicationNormalizeGovernors(host$cpu$governors)
    valid <- identical(preflight$schema_version, "1.1.0") &&
        identical(preflight$power_policy$mode, "observe_do_not_modify") &&
        identical(preflight$power_policy$mutation_allowed, FALSE) &&
        identical(
            governors,
            publicationNormalizeGovernors(preflight$power_policy$governors)
        ) && length(governors) > 0L &&
        publicationScalarNumber(
            maximum_run_allocated_disk_bytes,
            positive = TRUE
        )
    if (!valid) {
        publicationAbort(
            "Dynamic runtime safety authority is invalid",
            "multischolar_publication_host_preflight_error"
        )
    }
    list(
        safety_mode = "dynamic_observed_v1",
        minimum_available_memory_bytes = preflight$required_available_memory_bytes,
        minimum_available_disk_bytes = preflight$required_available_disk_bytes,
        maximum_load = host$cpu$logical_cores * 0.10,
        maximum_thermal_celsius = max(
            85,
            host$thermal$maximum_celsius + 10
        ),
        baseline_governors = as.list(governors),
        frequency_telemetry_required = TRUE,
        governor_stability_required = TRUE,
        maximum_run_allocated_disk_bytes = maximum_run_allocated_disk_bytes
    )
}

publicationValidateRuntimeSafetyLimits <- function(limits) {
    common <- c(
        "minimum_available_memory_bytes", "minimum_available_disk_bytes",
        "maximum_load", "maximum_thermal_celsius",
        "maximum_run_allocated_disk_bytes"
    )
    legacy <- c(
        common,
        "baseline_cpu_frequency_hz",
        "maximum_cpu_frequency_drift_fraction"
    )
    dynamic <- c(
        "safety_mode", common, "baseline_governors",
        "frequency_telemetry_required", "governor_stability_required"
    )
    mode <- if (is.list(limits) && setequal(names(limits), legacy)) {
        "legacy_frequency_drift_v1"
    } else if (is.list(limits) && setequal(names(limits), dynamic) &&
        identical(limits$safety_mode, "dynamic_observed_v1")) {
        "dynamic_observed_v1"
    } else {
        NULL
    }
    numeric_valid <- !is.null(mode) && all(vapply(
        limits[common],
        publicationScalarNumber,
        logical(1),
        positive = TRUE
    ))
    mode_valid <- if (identical(mode, "legacy_frequency_drift_v1")) {
        all(vapply(
            limits[setdiff(legacy, common)],
            publicationScalarNumber,
            logical(1),
            positive = TRUE
        ))
    } else if (identical(mode, "dynamic_observed_v1")) {
        length(publicationNormalizeGovernors(limits$baseline_governors)) > 0L &&
            isTRUE(limits$frequency_telemetry_required) &&
            isTRUE(limits$governor_stability_required)
    } else {
        FALSE
    }
    if (!numeric_valid || !mode_valid) {
        publicationAbort(
            "Runtime safety limits are invalid",
            "multischolar_publication_host_preflight_error"
        )
    }
    mode
}

publicationEvaluateRuntimeSafety <- function(sample, host, limits) {
    mode <- publicationValidateRuntimeSafetyLimits(limits)
    frequency <- host$frequencies$median_hz
    drift <- if (identical(mode, "legacy_frequency_drift_v1") &&
        publicationScalarNumber(frequency, positive = TRUE) &&
        limits$baseline_cpu_frequency_hz > 0) {
        abs(frequency - limits$baseline_cpu_frequency_hz) /
            limits$baseline_cpu_frequency_hz
    } else {
        Inf
    }
    load_values <- as.numeric(unlist(host$load_average, use.names = FALSE))
    load_valid <- length(load_values) == 3L && all(is.finite(load_values)) &&
        if (identical(mode, "dynamic_observed_v1")) {
            load_values[[1L]] <= limits$maximum_load
        } else {
            all(load_values <= limits$maximum_load)
        }
    common_checks <- list(
        cgroup_swap = publicationScalarNumber(sample$memory$swap_current_bytes) &&
            sample$memory$swap_current_bytes == 0,
        run_disk = publicationScalarNumber(sample$disk$allocated_bytes) &&
            sample$disk$allocated_bytes <= limits$maximum_run_allocated_disk_bytes,
        host_memory = publicationScalarNumber(host$memory$MemAvailable) &&
            host$memory$MemAvailable >= limits$minimum_available_memory_bytes,
        host_disk = isTRUE(host$filesystem$available) &&
            publicationScalarNumber(host$filesystem$available_bytes) &&
            host$filesystem$available_bytes >= limits$minimum_available_disk_bytes,
        host_load = load_valid,
        thermal = isTRUE(host$thermal$available) &&
            publicationScalarNumber(host$thermal$maximum_celsius) &&
            host$thermal$maximum_celsius <= limits$maximum_thermal_celsius
    )
    frequency_checks <- if (identical(mode, "dynamic_observed_v1")) {
        observed_governors <- publicationNormalizeGovernors(host$governors)
        frequency_valid <- isTRUE(host$frequencies$available) && all(vapply(
            host$frequencies[c("minimum_hz", "median_hz", "maximum_hz")],
            publicationScalarNumber,
            logical(1),
            positive = TRUE
        ))
        list(
            frequency_telemetry = frequency_valid,
            governor_stability = identical(
                observed_governors,
                publicationNormalizeGovernors(limits$baseline_governors)
            )
        )
    } else {
        list(
            frequency = publicationScalarNumber(drift) &&
                drift <= limits$maximum_cpu_frequency_drift_fraction
        )
    }
    checks <- c(common_checks, frequency_checks)
    safe <- isTRUE(all(unlist(checks, use.names = FALSE)))
    list(
        safe = safe,
        safety_mode = mode,
        checks = checks,
        frequency_drift_fraction = drift,
        observation = list(
            governors = as.list(publicationNormalizeGovernors(host$governors)),
            frequencies = host$frequencies,
            load_average = host$load_average,
            maximum_thermal_celsius = host$thermal$maximum_celsius,
            available_memory_bytes = host$memory$MemAvailable,
            available_disk_bytes = host$filesystem$available_bytes
        ),
        reason = if (safe) NULL else paste(
            names(checks)[!unlist(checks, use.names = FALSE)],
            collapse = ","
        )
    )
}

publicationSafetyTraceEvidence <- function(
    trace,
    timeout_seconds,
    interval_seconds
) {
    maximum_records <- ceiling(timeout_seconds / interval_seconds) + 1L
    valid <- is.list(trace) && length(trace) > 0L &&
        length(trace) <= maximum_records && all(vapply(trace, \(record) {
            is.list(record) && publicationScalarNumber(
                record$elapsed_seconds
            ) && is.logical(record$safe) && length(record$safe) == 1L &&
                is.list(record$checks) && is.list(record$observation)
        }, logical(1)))
    governors <- unique(unlist(lapply(trace, \(record) {
        publicationNormalizeGovernors(record$observation$governors)
    }), use.names = FALSE))
    medians <- as.numeric(unlist(lapply(trace, \(record) {
        record$observation$frequencies$median_hz
    }), use.names = FALSE))
    evidence <- list(
        schema = "multischolar.omics_publication_safety_trace",
        schema_version = "1.0.0",
        interval_seconds = interval_seconds,
        maximum_records = as.integer(maximum_records),
        record_count = length(trace),
        valid = valid,
        governor_set = as.list(sort(unique(governors), method = "radix")),
        median_frequency_hz = if (length(medians) && all(is.finite(medians))) {
            stats::median(medians)
        } else {
            NULL
        },
        trace = trace
    )
    basis <- evidence
    basis$trace_digest <- NULL
    evidence$trace_digest <- publicationObjectDigest(basis)
    evidence
}

publicationValidateSafetyTraceEvidence <- function(evidence) {
    expected <- c(
        "schema", "schema_version", "interval_seconds", "maximum_records",
        "record_count", "valid", "governor_set", "median_frequency_hz",
        "trace", "trace_digest"
    )
    publicationRequireNames(evidence, expected, "Publication safety trace")
    valid_header <- identical(
        evidence$schema,
        "multischolar.omics_publication_safety_trace"
    ) && identical(evidence$schema_version, "1.0.0") &&
        publicationScalarNumber(evidence$interval_seconds, positive = TRUE) &&
        publicationScalarNumber(evidence$maximum_records, positive = TRUE) &&
        identical(evidence$record_count, length(evidence$trace)) &&
        isTRUE(evidence$valid) && publicationScalarString(evidence$trace_digest) &&
        grepl("^[0-9a-f]{64}$", evidence$trace_digest)
    timeout_seconds <- (evidence$maximum_records - 1) *
        evidence$interval_seconds
    replay <- if (valid_header) {
        publicationSafetyTraceEvidence(
            evidence$trace,
            timeout_seconds,
            evidence$interval_seconds
        )
    } else {
        NULL
    }
    if (!valid_header || !identical(evidence, replay)) {
        publicationAbort(
            "Publication safety trace is invalid or digest-mismatched",
            "multischolar_publication_host_preflight_error"
        )
    }
    invisible(evidence)
}

publicationRuntimeSafetyMonitor <- function(
    limits,
    filesystem_path = .PUBLICATION_REPO_ROOT,
    interval_seconds = 1
) {
    mode <- publicationValidateRuntimeSafetyLimits(limits)
    last_checked <- -Inf
    last_result <- list(
        safe = TRUE,
        fresh = FALSE,
        checks = list(),
        reason = NULL
    )
    monitor <- function(sample, elapsed_seconds) {
        immediate <- list(
            safe = isTRUE(sample$memory$swap_current_bytes == 0) &&
                isTRUE(sample$disk$allocated_bytes <=
                    limits$maximum_run_allocated_disk_bytes),
            safety_mode = mode,
            checks = list(
                cgroup_swap = isTRUE(sample$memory$swap_current_bytes == 0),
                run_disk = isTRUE(sample$disk$allocated_bytes <=
                    limits$maximum_run_allocated_disk_bytes)
            ),
            observation = list(),
            reason = NULL
        )
        if (!immediate$safe) {
            immediate$reason <- "cgroup_swap_or_run_disk"
            immediate$fresh <- TRUE
            immediate$elapsed_seconds <- elapsed_seconds
            return(immediate)
        }
        if (elapsed_seconds - last_checked < interval_seconds) {
            cached <- last_result
            cached$fresh <- FALSE
            return(cached)
        }
        last_checked <<- elapsed_seconds
        last_result <<- publicationEvaluateRuntimeSafety(
            sample,
            publicationLightweightHostSnapshot(filesystem_path),
            limits
        )
        last_result$fresh <- TRUE
        last_result$elapsed_seconds <- elapsed_seconds
        last_result
    }
    attr(monitor, "publication_safety_interval_seconds") <- interval_seconds
    monitor
}

publicationPilotQualificationDecision <- function(
    measurement,
    memory_threshold_bytes = 4 * 1024^3,
    elapsed_threshold_seconds = 60
) {
    phase_valid <- is.list(measurement$phase_evidence) &&
        isTRUE(measurement$phase_evidence$valid)
    safety_valid <- tryCatch({
        publicationValidateSafetyTraceEvidence(measurement$safety_evidence)
        TRUE
    }, error = \(...) FALSE)
    decision_valid <- identical(measurement$status, "passed") &&
        isTRUE(measurement$publication_certifiable) &&
        !isTRUE(measurement$timed_out) &&
        !isTRUE(measurement$safety_aborted) && phase_valid && safety_valid &&
        publicationScalarNumber(
            measurement$metrics$peak_charged_memory_bytes
        ) && publicationScalarNumber(measurement$metrics$elapsed_seconds)
    qualified <- decision_valid &&
        (measurement$metrics$peak_charged_memory_bytes >=
            memory_threshold_bytes ||
            measurement$metrics$elapsed_seconds >= elapsed_threshold_seconds)
    list(
        decision_valid = decision_valid,
        dimension_selection_allowed = decision_valid,
        qualified = qualified,
        decision = if (!decision_valid) {
            "none"
        } else if (qualified) {
            "qualified"
        } else {
            "not_qualified"
        },
        memory_threshold_bytes = memory_threshold_bytes,
        elapsed_threshold_seconds = elapsed_threshold_seconds,
        observed_peak_charged_memory_bytes =
            measurement$metrics$peak_charged_memory_bytes,
        observed_elapsed_seconds = measurement$metrics$elapsed_seconds,
        criterion = "peak_charged_memory_at_least_4_gib_or_elapsed_at_least_60_s"
    )
}

publicationPilotTerminalStatus <- function(measurement, qualification) {
    if (isTRUE(qualification$decision_valid)) {
        if (isTRUE(qualification$qualified)) return("qualified")
        return("not_qualified")
    }
    if (isTRUE(measurement$safety_aborted)) {
        return("safety_aborted_no_dimension_decision")
    }
    if (isTRUE(measurement$timed_out)) {
        return("timed_out_no_dimension_decision")
    }
    if (is.list(measurement$phase_evidence) &&
        identical(measurement$phase_evidence$valid, FALSE)) {
        return("phase_invalid_no_dimension_decision")
    }
    if (!is.list(measurement$safety_evidence) ||
        !isTRUE(measurement$safety_evidence$valid)) {
        return("safety_evidence_invalid_no_dimension_decision")
    }
    "worker_failed_no_dimension_decision"
}
