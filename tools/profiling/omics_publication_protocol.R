publicationScriptRepoRoot <- function(start = getwd()) {
    path <- normalizePath(start, mustWork = TRUE)
    repeat {
        if (file.exists(file.path(path, "DESCRIPTION")) &&
            dir.exists(file.path(path, "tools", "profiling"))) {
            return(path)
        }
        parent <- dirname(path)
        if (identical(parent, path)) {
            stop("Cannot locate MultiScholaR repository root", call. = FALSE)
        }
        path <- parent
    }
}

.PUBLICATION_REPO_ROOT <- publicationScriptRepoRoot()

publicationAbort <- function(message, class) {
    rlang::abort(
        message,
        class = c(class, "multischolar_publication_error")
    )
}

publicationUtcNow <- function() {
    format(
        as.POSIXct(Sys.time(), tz = "UTC"),
        "%Y-%m-%dT%H:%M:%OS3Z",
        tz = "UTC"
    )
}

publicationPath <- function(...) {
    file.path(.PUBLICATION_REPO_ROOT, ...)
}

publicationReadJson <- function(path) {
    resolved <- if (startsWith(path, .Platform$file.sep)) {
        path
    } else {
        publicationPath(path)
    }
    if (!file.exists(resolved)) {
        publicationAbort(
            paste("Publication JSON does not exist:", path),
            "multischolar_publication_missing_record"
        )
    }
    jsonlite::read_json(resolved, simplifyVector = FALSE)
}

publicationWriteJson <- function(value, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    temporary <- tempfile("publication-json-", tmpdir = dirname(path))
    on.exit(unlink(temporary, force = TRUE), add = TRUE)
    jsonlite::write_json(
        value,
        temporary,
        auto_unbox = TRUE,
        pretty = TRUE,
        null = "null",
        na = "null",
        digits = 17
    )
    if (!file.rename(temporary, path)) {
        publicationAbort(
            paste("Could not atomically publish JSON:", path),
            "multischolar_publication_write_failed"
        )
    }
    invisible(path)
}

publicationFileDigest <- function(path) {
    resolved <- if (startsWith(path, .Platform$file.sep)) {
        path
    } else {
        publicationPath(path)
    }
    digest::digest(file = resolved, algo = "sha256", serialize = FALSE)
}

publicationSortObject <- function(value) {
    if (!is.list(value)) return(value)
    if (!is.null(names(value))) {
        value <- value[order(names(value), method = "radix")]
    }
    lapply(value, publicationSortObject)
}

publicationObjectDigest <- function(value) {
    canonical <- publicationSortObject(value)
    encoded <- jsonlite::toJSON(
        canonical,
        auto_unbox = TRUE,
        null = "null",
        na = "null",
        digits = 17,
        pretty = FALSE
    )
    digest::digest(encoded, algo = "sha256", serialize = FALSE)
}

publicationWithPreservedSeed <- function(seed, code) {
    seed_existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    previous_seed <- if (seed_existed) get(".Random.seed", envir = .GlobalEnv) else
        NULL
    on.exit({
        if (seed_existed) {
            assign(".Random.seed", previous_seed, envir = .GlobalEnv)
        } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            rm(".Random.seed", envir = .GlobalEnv)
        }
    }, add = TRUE)
    set.seed(seed)
    force(code)
}

publicationMetricContractPath <- function() {
    "tests/testdata/omics-performance/metric-contract-v1.json"
}

publicationValidateMetricContract <- function(contract) {
    publicationRequireNames(contract, c(
        "schema", "schema_version", "metric_contract_id", "owner_ticket_id",
        "measurement_schema", "authoritative_references", "rules", "metrics"
    ), "Metric contract")
    required_metrics <- c(
        "peak_charged_memory_bytes", "retained_charged_memory_bytes",
        "peak_anonymous_memory_bytes", "retained_anonymous_memory_bytes",
        "peak_file_memory_bytes", "retained_file_memory_bytes",
        "peak_kernel_memory_bytes", "retained_kernel_memory_bytes",
        "peak_slab_memory_bytes", "retained_slab_memory_bytes",
        "peak_swap_bytes", "peak_pss_bytes", "retained_pss_bytes",
        "peak_uss_bytes", "retained_uss_bytes", "peak_rss_bytes",
        "retained_rss_bytes", "elapsed_seconds", "phase_cpu_seconds",
        "cpu_usage_seconds", "primary_work_units_per_wall_second",
        "primary_work_units_per_cpu_second", "io", "phase_io",
        "peak_logical_disk_bytes", "peak_allocated_disk_bytes",
        "final_logical_disk_bytes", "final_allocated_disk_bytes",
        "final_file_count", "query_p95_seconds", "arrow_pool_bytes",
        "duckdb_resources"
    )
    metric_fields <- c(
        "metric_id", "unit", "source", "aggregation", "phase", "role"
    )
    metric_ids <- vapply(contract$metrics, \(metric) {
        publicationRequireNames(metric, metric_fields, "Metric definition")
        fields_valid <- all(vapply(metric, publicationScalarString, logical(1)))
        if (!fields_valid) {
            publicationAbort(
                "Metric definition contains an empty field",
                "multischolar_publication_schema_error"
            )
        }
        metric$metric_id
    }, character(1))
    references <- unlist(contract$authoritative_references, use.names = FALSE)
    valid <- identical(
        contract$schema,
        "multischolar.omics_publication_metric_contract"
    ) && identical(contract$schema_version, "1.0.0") &&
        identical(contract$owner_ticket_id, "OMICS-ART-063") &&
        identical(
            contract$measurement_schema,
            "multischolar.omics_publication_cgroup_measurement/1.0.0"
        ) && !anyDuplicated(metric_ids) &&
        all(required_metrics %in% metric_ids) && length(references) >= 3L &&
        all(startsWith(references, "https://www.kernel.org/")) &&
        identical(contract$rules$absent_is_zero, FALSE) &&
        identical(contract$rules$measurement_overhead_subtracted, FALSE) &&
        identical(contract$rules$memory_stat_slab_is_subset_of_kernel, TRUE) &&
        identical(contract$rules$retention_acknowledgement, "fifo_v1") &&
        identical(contract$rules$maximum_boundary_bracket_seconds, 0.5) &&
        identical(contract$rules$query_resampling_unit, "fresh_process_run_pair")
    if (!valid) {
        publicationAbort(
            "Metric contract differs from publication semantics",
            "multischolar_publication_schema_error"
        )
    }
    invisible(contract)
}

publicationMetricContractBinding <- function() {
    path <- publicationMetricContractPath()
    contract <- publicationReadJson(path)
    publicationValidateMetricContract(contract)
    list(
        metric_contract_id = contract$metric_contract_id,
        path = path,
        sha256 = publicationFileDigest(path)
    )
}

publicationKernelSourcePaths <- function() {
    c(
        "tools/profiling/omics_publication_protocol.R",
        "tools/profiling/omics_publication_measure_spec.R",
        "tools/profiling/omics_publication_linux_resources.R",
        "tools/profiling/omics_publication_retained_resources.R",
        "tools/profiling/omics_publication_host_safety.R",
        "tools/profiling/omics_publication_schedule.R",
        "tools/profiling/omics_publication_campaign_state.R",
        "tools/profiling/omics_publication_statistics.R",
        "tools/profiling/run_omics_publication_benchmark.R",
        "tools/profiling/omics_publication_oom_probe.c"
    )
}

publicationKernelSourceBindings <- function() {
    lapply(publicationKernelSourcePaths(), \(path) {
        list(path = path, sha256 = publicationFileDigest(path))
    })
}

publicationDynamicFrequencySuccessorPath <- function() {
    "tests/testdata/omics-performance/dynamic-frequency-successor-v1.json"
}

publicationBindingMatchesCurrentOrSuccessor <- function(binding) {
    resolved_path <- if (is.list(binding) &&
        publicationScalarString(binding$path) &&
        startsWith(binding$path, .Platform$file.sep)) {
        binding$path
    } else if (is.list(binding) && publicationScalarString(binding$path)) {
        publicationPath(binding$path)
    } else {
        NULL
    }
    if (!is.list(binding) ||
        !setequal(names(binding), c("path", "sha256")) ||
        !publicationScalarString(binding$path) ||
        !publicationScalarString(binding$sha256) ||
        !grepl("^[0-9a-f]{64}$", binding$sha256) ||
        is.null(resolved_path) || !file.exists(resolved_path)) {
        return(FALSE)
    }
    current_sha256 <- publicationFileDigest(binding$path)
    if (identical(current_sha256, binding$sha256)) return(TRUE)
    successor_path <- publicationDynamicFrequencySuccessorPath()
    if (!file.exists(publicationPath(successor_path))) return(FALSE)
    successor <- publicationReadJson(successor_path)
    predecessor_valid <- identical(
        publicationFileDigest(successor$predecessor$path),
        successor$predecessor$sha256
    )
    successor_valid <- identical(
        publicationFileDigest(successor$successor$path),
        successor$successor$sha256
    )
    transitions <- successor$kernel_transitions
    matches <- Filter(\(transition) {
        identical(transition$path, binding$path) &&
            identical(transition$predecessor_sha256, binding$sha256) &&
            identical(transition$current_sha256, current_sha256)
    }, transitions)
    isTRUE(predecessor_valid) && isTRUE(successor_valid) &&
        length(matches) == 1L
}

publicationValidateKernelSourceBindings <- function(bindings) {
    expected <- publicationKernelSourcePaths()
    paths <- vapply(bindings, `[[`, character(1), "path")
    if (!identical(paths, expected) || any(!vapply(
        bindings,
        publicationBindingMatchesCurrentOrSuccessor,
        logical(1)
    ))) {
        publicationAbort(
            "Publication kernel source binding differs",
            "multischolar_publication_binding_error"
        )
    }
    invisible(bindings)
}

publicationScalarString <- function(value) {
    is.character(value) && length(value) == 1L && !is.na(value) && nzchar(value)
}

publicationScalarNumber <- function(value, positive = FALSE) {
    valid <- is.numeric(value) && length(value) == 1L && is.finite(value)
    valid && (!positive || value > 0)
}

publicationRequireNames <- function(value, expected, label) {
    if (!is.list(value) || !setequal(names(value), expected)) {
        publicationAbort(
            paste(label, "fields differ from the frozen contract"),
            "multischolar_publication_schema_error"
        )
    }
    invisible(value)
}

publicationReadGovernanceManifest <- function(
    path = "tests/testdata/omics-performance/governance-manifest-v1.json"
) {
    manifest <- publicationReadJson(path)
    expected <- c(
        "schema", "schema_version", "manifest_id", "owner_ticket_id",
        "records", "immutability"
    )
    publicationRequireNames(manifest, expected, "Governance manifest")
    if (!identical(
        manifest$schema,
        "multischolar.omics_publication_governance_manifest"
    ) || !identical(manifest$schema_version, "1.0.0")) {
        publicationAbort(
            "Governance manifest version is unsupported",
            "multischolar_publication_schema_error"
        )
    }
    paths <- vapply(manifest$records, `[[`, character(1), "path")
    if (anyDuplicated(paths) || any(!nzchar(paths))) {
        publicationAbort(
            "Governance manifest contains duplicate record paths",
            "multischolar_publication_binding_error"
        )
    }
    for (record in manifest$records) {
        if (!identical(publicationFileDigest(record$path), record$sha256)) {
            publicationAbort(
                paste("Governance record digest mismatch:", record$path),
                "multischolar_publication_binding_error"
            )
        }
    }
    manifest
}

publicationValidateGovernanceManifestV2 <- function(record) {
    expected <- c(
        "schema", "schema_version", "manifest_id", "owner_ticket_id",
        "predecessor", "records", "implementation_sources", "decisions",
        "immutability"
    )
    publicationRequireNames(record, expected, "Governance manifest v2")
    if (!identical(
        record$schema,
        "multischolar.omics_publication_governance_manifest"
    ) || !identical(record$schema_version, "2.0.0") ||
        !identical(record$owner_ticket_id, "OMICS-ART-063")) {
        publicationAbort(
            "Governance manifest v2 identity differs",
            "multischolar_publication_schema_error"
        )
    }
    bindings <- c(list(record$predecessor), record$records, record$implementation_sources)
    paths <- vapply(bindings, `[[`, character(1), "path")
    if (anyDuplicated(paths)) {
        publicationAbort(
            "Governance manifest v2 paths are duplicate",
            "multischolar_publication_binding_error"
        )
    }
    for (binding in bindings) {
        if (!identical(publicationFileDigest(binding$path), binding$sha256)) {
            publicationAbort(
                paste("Governance manifest v2 digest mismatch:", binding$path),
                "multischolar_publication_binding_error"
            )
        }
    }
    decisions <- record$decisions
    if (!identical(decisions$exact_primary_pair_count, 48L) ||
        !identical(decisions$pair_count_status, "precision_satisfied") ||
        isTRUE(decisions$candidate_loaded_for_precision) ||
        isTRUE(decisions$publication_host_certified) ||
        isTRUE(decisions$campaign_execution_authorized) ||
        isTRUE(decisions$promotion_authorized) ||
        !isTRUE(record$immutability$additive_successor_required) ||
        isTRUE(record$immutability$mutate_version_in_place)) {
        publicationAbort(
            "Governance manifest v2 decisions are unsafe",
            "multischolar_publication_binding_error"
        )
    }
    invisible(record)
}

publicationValidateGovernanceManifestV3 <- function(record) {
    expected <- c(
        "schema", "schema_version", "manifest_id", "owner_ticket_id",
        "predecessor", "records", "implementation_sources", "decisions",
        "immutability"
    )
    publicationRequireNames(record, expected, "Governance manifest v3")
    identity_valid <- identical(
        record$schema,
        "multischolar.omics_publication_governance_manifest"
    ) && identical(record$schema_version, "3.0.0") &&
        identical(record$owner_ticket_id, "OMICS-ART-063") &&
        identical(
            record$predecessor$path,
            "tests/testdata/omics-performance/governance-manifest-v2.json"
        )
    if (!identity_valid) {
        publicationAbort(
            "Governance manifest v3 identity differs",
            "multischolar_publication_schema_error"
        )
    }
    immutable_bindings <- c(list(record$predecessor), record$records)
    immutable_paths <- vapply(
        immutable_bindings,
        `[[`,
        character(1),
        "path"
    )
    immutable_valid <- !anyDuplicated(immutable_paths) && all(vapply(
        immutable_bindings,
        \(binding) {
        publicationScalarString(binding$sha256) &&
            grepl("^[0-9a-f]{64}$", binding$sha256) &&
            identical(publicationFileDigest(binding$path), binding$sha256)
        },
        logical(1)
    ))
    record_paths <- stats::setNames(
        vapply(record$records, `[[`, character(1), "path"),
        vapply(record$records, `[[`, character(1), "record_id")
    )
    expected_record_ids <- c(
        "host_preflight_contract", "metric_contract", "cache_contract",
        "null_cgroup_calibration", "precision_successor"
    )
    source_paths <- vapply(
        record$implementation_sources,
        `[[`,
        character(1),
        "path"
    )
    sources_valid <- identical(source_paths, publicationKernelSourcePaths()) &&
        all(vapply(
            record$implementation_sources,
            publicationBindingMatchesCurrentOrSuccessor,
            logical(1)
        ))
    if (!immutable_valid ||
        !identical(names(record_paths), expected_record_ids) ||
        !sources_valid) {
        publicationAbort(
            "Governance manifest v3 binding differs",
            "multischolar_publication_binding_error"
        )
    }
    publicationValidateHostPreflightContract(publicationReadJson(
        record_paths[["host_preflight_contract"]]
    ))
    publicationValidateMetricContract(publicationReadJson(
        record_paths[["metric_contract"]]
    ))
    publicationValidateCacheContract(publicationReadJson(
        record_paths[["cache_contract"]]
    ))
    null_record <- publicationReadJson(record_paths[["null_cgroup_calibration"]])
    precision <- publicationReadJson(record_paths[["precision_successor"]])
    publicationValidateNullCalibration(null_record)
    publicationValidatePrecisionSuccessor(precision, null_record)
    decisions <- record$decisions
    decisions_valid <- identical(
        decisions$exact_primary_pair_count,
        precision$selected_pairs
    ) && identical(decisions$pair_count_status, precision$status) &&
        !isTRUE(decisions$candidate_loaded_for_precision) &&
        !isTRUE(decisions$publication_host_certified) &&
        !isTRUE(decisions$campaign_execution_authorized) &&
        !isTRUE(decisions$promotion_authorized) &&
        isTRUE(record$immutability$additive_successor_required) &&
        !isTRUE(record$immutability$mutate_version_in_place)
    if (!decisions_valid) {
        publicationAbort(
            "Governance manifest v3 decisions are unsafe",
            "multischolar_publication_binding_error"
        )
    }
    invisible(record)
}

publicationGovernanceRecord <- function(record_id, manifest = NULL) {
    if (is.null(manifest)) manifest <- publicationReadGovernanceManifest()
    ids <- vapply(manifest$records, `[[`, character(1), "record_id")
    matches <- which(ids == record_id)
    if (length(matches) != 1L) {
        publicationAbort(
            paste("Governance record is not uniquely owned:", record_id),
            "multischolar_publication_binding_error"
        )
    }
    publicationReadJson(manifest$records[[matches]]$path)
}

publicationProtocolBundle <- function() {
    manifest <- publicationReadGovernanceManifest()
    list(
        manifest = manifest,
        protocol = publicationGovernanceRecord("protocol", manifest),
        coverage = publicationGovernanceRecord("coverage", manifest),
        roles = publicationGovernanceRecord("roles", manifest),
        projects = publicationGovernanceRecord("projects", manifest),
        estimands = publicationGovernanceRecord("estimands", manifest),
        splits = publicationGovernanceRecord("splits", manifest),
        threshold_grid = publicationGovernanceRecord("threshold_grid", manifest),
        retry_policy = publicationGovernanceRecord("retry_policy", manifest),
        campaign_matrix = publicationGovernanceRecord("campaign_matrix", manifest),
        campaign_budget = publicationGovernanceRecord("campaign_budget", manifest),
        receipt_schema = publicationGovernanceRecord(
            "policy_receipt_schema",
            manifest
        )
    )
}

publicationGitRevision <- function() {
    result <- processx::run(
        "git",
        c("rev-parse", "HEAD"),
        wd = .PUBLICATION_REPO_ROOT,
        error_on_status = FALSE,
        echo = FALSE
    )
    if (result$status != 0L) return(NULL)
    trimws(result$stdout)
}

publicationPackageVersions <- function(packages) {
    values <- lapply(packages, \(package) {
        if (!requireNamespace(package, quietly = TRUE)) return(NULL)
        as.character(utils::packageVersion(package))
    })
    stats::setNames(values, packages)
}
