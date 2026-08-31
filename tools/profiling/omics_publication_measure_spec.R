publicationCacheContractPath <- function() {
    "tests/testdata/omics-performance/cache-contract-v1.json"
}

publicationValidateCacheContract <- function(contract) {
    publicationRequireNames(contract, c(
        "schema", "schema_version", "cache_contract_id", "owner_ticket_id",
        "strata", "rules"
    ), "Cache contract")
    strata_fields <- c(
        "stratum", "pre_read_required", "verified_page_cache_drop_required",
        "primary_comparison_eligible"
    )
    strata <- vapply(contract$strata, \(stratum) {
        publicationRequireNames(stratum, strata_fields, "Cache stratum")
        stratum$stratum
    }, character(1))
    valid <- identical(
        contract$schema,
        "multischolar.omics_publication_cache_contract"
    ) && identical(contract$schema_version, "1.0.0") &&
        identical(contract$owner_ticket_id, "OMICS-ART-063") &&
        setequal(strata, c(
            "standardized_warm_input", "first_open", "cold_uncontrolled",
            "verified_cold"
        )) && !anyDuplicated(strata) &&
        identical(contract$rules$bare_cold_label_allowed, FALSE) &&
        identical(
            contract$rules$cold_requires_privileged_drop_receipt,
            TRUE
        ) && identical(contract$rules$cache_strata_pooled, FALSE) &&
        identical(
            contract$rules$cache_evidence_alone_authorizes_promotion,
            FALSE
        )
    if (!valid) {
        publicationAbort(
            "Cache contract differs from publication semantics",
            "multischolar_publication_schema_error"
        )
    }
    invisible(contract)
}

publicationValidateCacheEvidence <- function(evidence, contract) {
    publicationValidateCacheContract(contract)
    publicationRequireNames(evidence, c(
        "cache_contract_id", "stratum", "input_sha256", "pre_read_complete",
        "page_cache_drop_attempted", "page_cache_drop_verified",
        "page_cache_drop_method", "page_cache_drop_receipt_sha256",
        "primary_comparison_eligible"
    ), "Cache evidence")
    if (!publicationScalarString(evidence$stratum)) {
        publicationAbort(
            "Cache stratum is invalid",
            "multischolar_publication_cache_error"
        )
    }
    strata <- stats::setNames(contract$strata, vapply(
        contract$strata,
        `[[`,
        character(1),
        "stratum"
    ))
    selected <- strata[[evidence$stratum]]
    digest_valid <- publicationScalarString(evidence$input_sha256) &&
        grepl("^[0-9a-f]{64}$", evidence$input_sha256)
    common_valid <- !is.null(selected) && digest_valid &&
        identical(evidence$cache_contract_id, contract$cache_contract_id) &&
        identical(
            evidence$primary_comparison_eligible,
            selected$primary_comparison_eligible
        )
    if (!common_valid || identical(evidence$stratum, "cold")) {
        publicationAbort(
            "Cache evidence binding is invalid",
            "multischolar_publication_cache_error"
        )
    }
    if (isTRUE(selected$pre_read_required) &&
        !isTRUE(evidence$pre_read_complete)) {
        publicationAbort(
            "Warm-input cache evidence lacks pre-read completion",
            "multischolar_publication_cache_error"
        )
    }
    if (isTRUE(selected$verified_page_cache_drop_required)) {
        receipt_valid <- publicationScalarString(
            evidence$page_cache_drop_receipt_sha256
        ) && grepl(
            "^[0-9a-f]{64}$",
            evidence$page_cache_drop_receipt_sha256
        )
        if (!isTRUE(evidence$page_cache_drop_attempted) ||
            !isTRUE(evidence$page_cache_drop_verified) ||
            !publicationScalarString(evidence$page_cache_drop_method) ||
            !receipt_valid) {
            publicationAbort(
                "Verified-cold cache evidence lacks a privileged receipt",
                "multischolar_publication_cache_error"
            )
        }
    } else if (isTRUE(evidence$page_cache_drop_verified) ||
        !is.null(evidence$page_cache_drop_receipt_sha256)) {
        publicationAbort(
            "Non-cold cache evidence contains a false drop claim",
            "multischolar_publication_cache_error"
        )
    }
    invisible(evidence)
}

publicationValidateSha256 <- function(value, label) {
    if (!publicationScalarString(value) || !grepl("^[0-9a-f]{64}$", value)) {
        publicationAbort(
            paste(label, "must be a lowercase SHA-256 digest"),
            "multischolar_publication_binding_error"
        )
    }
    invisible(value)
}

publicationValidateMeasureSpecHeader <- function(spec) {
    fields <- c(
        "schema", "schema_version", "measure_spec_id", "protocol_binding",
        "metric_contract", "estimand_binding", "schedule_binding",
        "source_binding", "candidate_binding", "cache_evidence", "work",
        "command", "arguments", "run_dir", "execution", "environment",
        "host_preflight", "host_preflight_sha256", "safety_limits"
    )
    publicationRequireNames(spec, fields, "Publication measure spec")
    valid <- identical(
        spec$schema,
        "multischolar.omics_publication_measure_spec"
    ) && identical(spec$schema_version, "1.0.0") &&
        publicationScalarString(spec$measure_spec_id) &&
        publicationScalarString(spec$command) &&
        publicationScalarString(spec$run_dir) &&
        startsWith(spec$run_dir, .Platform$file.sep) &&
        is.list(spec$arguments) &&
        all(vapply(spec$arguments, publicationScalarString, logical(1))) &&
        is.list(spec$execution) && is.list(spec$environment) &&
        isTRUE(spec$host_preflight$certified)
    if (!valid) {
        publicationAbort(
            "Publication measure spec fields are invalid",
            "multischolar_publication_measure_spec_error"
        )
    }
    invisible(spec)
}

publicationValidateMeasureGovernanceBindings <- function(spec) {
    protocol_path <- "tests/testdata/omics-performance/protocol-v1.json"
    protocol <- publicationReadJson(protocol_path)
    publicationRequireNames(spec$protocol_binding, c(
        "protocol_id", "path", "sha256"
    ), "Protocol binding")
    protocol_valid <- identical(
        spec$protocol_binding$protocol_id,
        protocol$protocol_id
    ) && identical(spec$protocol_binding$path, protocol_path) &&
        identical(
            spec$protocol_binding$sha256,
            publicationFileDigest(protocol_path)
        )
    if (!protocol_valid ||
        !identical(spec$metric_contract, publicationMetricContractBinding())) {
        publicationAbort(
            "Measure spec governance binding differs",
            "multischolar_publication_binding_error"
        )
    }
    cache_contract <- publicationReadJson(publicationCacheContractPath())
    publicationValidateCacheEvidence(spec$cache_evidence, cache_contract)
    invisible(spec)
}

publicationMeasureEntityBindingFields <- function() {
    list(
        estimand_binding = c("estimand_id", "estimands_sha256"),
        schedule_binding = c(
            "schedule_id", "schedule_sha256", "slot_id", "slot_sha256",
            "work_binding_sha256"
        ),
        source_binding = c(
            "project_id", "source_id", "source_sha256", "input_sha256",
            "exact_input_bytes", "evidence_class", "private_source",
            "private_values_retained", "privacy_receipt_sha256"
        ),
        candidate_binding = c("revision", "package_sha256")
    )
}

publicationValidateMeasureEntityBindings <- function(spec) {
    group_fields <- publicationMeasureEntityBindingFields()
    for (group in names(group_fields)) {
        publicationRequireNames(
            spec[[group]],
            group_fields[[group]],
            paste("Measure spec", group)
        )
    }
    identity_fields <- list(
        estimand_binding = "estimand_id",
        schedule_binding = c("schedule_id", "slot_id"),
        source_binding = c("project_id", "source_id")
    )
    identities_valid <- all(vapply(names(identity_fields), \(group) {
        all(vapply(
            spec[[group]][identity_fields[[group]]],
            publicationScalarString,
            logical(1)
        ))
    }, logical(1)))
    revision_valid <- publicationScalarString(spec$candidate_binding$revision) &&
        grepl("^[0-9a-f]{40}$", spec$candidate_binding$revision)
    if (!identities_valid || !revision_valid) {
        publicationAbort(
            "Measure spec entity identity is invalid",
            "multischolar_publication_binding_error"
        )
    }
    invisible(spec)
}

publicationValidateMeasurePrivacyBinding <- function(source) {
    valid <- source$evidence_class %in% c(
        "public_real", "private_calibrated", "public_generated",
        "fixture_correctness"
    ) && is.logical(source$private_source) &&
        length(source$private_source) == 1L &&
        is.logical(source$private_values_retained) &&
        length(source$private_values_retained) == 1L &&
        !isTRUE(source$private_values_retained) &&
        (!isTRUE(source$private_source) || identical(
            source$evidence_class,
            "private_calibrated"
        ))
    if (!valid) {
        publicationAbort(
            "Measure spec privacy binding is invalid",
            "multischolar_publication_binding_error"
        )
    }
    invisible(source)
}

publicationValidateMeasureBindingDigests <- function(spec) {
    groups <- names(publicationMeasureEntityBindingFields())
    for (group in groups) {
        digests <- names(spec[[group]])[grepl("sha256$", names(spec[[group]]))]
        if (!length(digests)) {
            publicationAbort(
                paste("Measure spec", group, "lacks a digest"),
                "multischolar_publication_binding_error"
            )
        }
        for (name in digests) {
            publicationValidateSha256(spec[[group]][[name]], paste(group, name))
        }
    }
    publicationValidateSha256(
        spec$host_preflight_sha256,
        "host_preflight_sha256"
    )
    host_valid <- identical(
        spec$host_preflight_sha256,
        publicationObjectDigest(spec$host_preflight)
    ) && identical(
        spec$host_preflight$schema,
        "multischolar.omics_publication_host_preflight"
    )
    if (!host_valid) {
        publicationAbort(
            "Measure spec host preflight binding differs",
            "multischolar_publication_binding_error"
        )
    }
    invisible(spec)
}

publicationValidateMeasureWorkBinding <- function(spec) {
    fields <- c(
        "primary_work_unit_id", "expected_primary_work_count",
        "exact_input_bytes"
    )
    publicationRequireNames(spec$work, fields, "Measure work binding")
    work_valid <- publicationScalarString(spec$work$primary_work_unit_id) &&
        publicationScalarNumber(
            spec$work$expected_primary_work_count,
            positive = TRUE
        ) && publicationScalarNumber(spec$work$exact_input_bytes, positive = TRUE)
    bindings_valid <- identical(
        spec$schedule_binding$work_binding_sha256,
        publicationObjectDigest(spec$work)
    ) && identical(
        spec$cache_evidence$input_sha256,
        spec$source_binding$input_sha256
    ) && identical(
        spec$work$exact_input_bytes,
        spec$source_binding$exact_input_bytes
    )
    if (!work_valid || !bindings_valid) {
        publicationAbort(
            "Measure spec source cache schedule and work bindings differ",
            "multischolar_publication_binding_error"
        )
    }
    invisible(spec)
}

publicationValidateMeasureRuntime <- function(spec) {
    thread_names <- c(
        "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
        "ARROW_NUM_THREADS", "DUCKDB_THREADS"
    )
    environment <- unlist(spec$environment, use.names = TRUE)
    execution_fields <- c(
        "sampling_interval_ms", "disk_sampling_interval_ms", "timeout_seconds",
        "retained_dwell_seconds", "retained_window_seconds",
        "retention_acknowledgement", "maximum_boundary_bracket_seconds",
        "retained_boundary_tolerance_ms"
    )
    publicationRequireNames(spec$execution, execution_fields, "Measure execution")
    execution <- unlist(spec$execution, use.names = TRUE)
    numeric_execution <- execution[setdiff(
        names(execution),
        "retention_acknowledgement"
    )]
    execution_valid <- all(is.finite(as.numeric(numeric_execution))) &&
        execution[["sampling_interval_ms"]] == 20 &&
        execution[["disk_sampling_interval_ms"]] == 500 &&
        execution[["retained_dwell_seconds"]] == 5 &&
        execution[["retained_window_seconds"]] == 2 &&
        execution[["retention_acknowledgement"]] == "fifo_v1" &&
        execution[["maximum_boundary_bracket_seconds"]] == 0.5 &&
        execution[["retained_boundary_tolerance_ms"]] == 100 &&
        execution[["timeout_seconds"]] > 5
    threads_valid <- all(thread_names %in% names(environment)) &&
        all(environment[thread_names] == "1")
    safety_mode <- tryCatch({
        publicationValidateRuntimeSafetyLimits(spec$safety_limits)
    }, error = \(...) NULL)
    safety_valid <- !is.null(safety_mode)
    if (identical(safety_mode, "dynamic_observed_v1")) {
        safety_valid <- identical(
            spec$host_preflight$schema_version,
            "1.1.0"
        ) && identical(
            spec$host_preflight$power_policy$mode,
            "observe_do_not_modify"
        ) && identical(
            spec$host_preflight$power_policy$mutation_allowed,
            FALSE
        ) && identical(
            publicationNormalizeGovernors(
                spec$safety_limits$baseline_governors
            ),
            publicationNormalizeGovernors(
                spec$host_preflight$power_policy$governors
            )
        )
    }
    if (!execution_valid || !threads_valid ||
        !safety_valid) {
        publicationAbort(
            "Measure runtime differs from the publication protocol",
            "multischolar_publication_measure_spec_error"
        )
    }
    invisible(spec)
}

publicationValidateMeasureSpec <- function(spec) {
    publicationValidateMeasureSpecHeader(spec)
    publicationValidateMeasureGovernanceBindings(spec)
    publicationValidateMeasureEntityBindings(spec)
    publicationValidateMeasurePrivacyBinding(spec$source_binding)
    publicationValidateMeasureBindingDigests(spec)
    publicationValidateMeasureWorkBinding(spec)
    publicationValidateMeasureRuntime(spec)
    invisible(spec)
}
