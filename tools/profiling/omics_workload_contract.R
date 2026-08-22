.OMICS_WORKLOAD_SCHEMA <- "multischolar.omics_workload_contract"
.OMICS_WORKLOAD_SCHEMA_VERSION <- "1.0.0"
.OMICS_WORKLOAD_EVIDENCE_SCHEMA <- "multischolar.omics_workload_evidence"
.OMICS_WORKLOAD_EVIDENCE_VERSION <- "1.0.0"

omicsWorkloadAbort <- function(message, class = "omics_workload_contract_error") {
    condition <- structure(
        list(message = message, call = NULL),
        class = c(class, "error", "condition")
    )
    stop(condition)
}

omicsWorkloadRequireNames <- function(value, expected, label) {
    if (!is.list(value) || is.null(names(value))) {
        omicsWorkloadAbort(sprintf("%s must be a named list", label))
    }
    missing <- setdiff(expected, names(value))
    unknown <- setdiff(names(value), expected)
    if (length(missing) || length(unknown)) {
        details <- c(
            if (length(missing)) sprintf("missing: %s", paste(missing, collapse = ", ")),
            if (length(unknown)) sprintf("unknown: %s", paste(unknown, collapse = ", "))
        )
        omicsWorkloadAbort(sprintf(
            "%s fields are invalid (%s)",
            label,
            paste(details, collapse = "; ")
        ))
    }
    invisible(value)
}

omicsWorkloadRequireText <- function(value, label) {
    if (!is.character(value) || length(value) != 1L ||
        is.na(value) || !nzchar(value)) {
        omicsWorkloadAbort(sprintf("%s must be one non-empty string", label))
    }
    invisible(value)
}

omicsWorkloadRequireCount <- function(value, label, minimum = 1L) {
    valid <- is.numeric(value) && length(value) == 1L && is.finite(value) &&
        value == as.integer(value) && value >= minimum
    if (!valid) {
        omicsWorkloadAbort(sprintf(
            "%s must be one integer greater than or equal to %d",
            label,
            minimum
        ))
    }
    invisible(value)
}

omicsWorkloadRequireDigest <- function(value, label) {
    omicsWorkloadRequireText(value, label)
    if (!grepl("^[a-f0-9]{64}$", value)) {
        omicsWorkloadAbort(sprintf("%s must be a lowercase SHA-256 digest", label))
    }
    invisible(value)
}

omicsWorkloadForbiddenPrivateFields <- function(value, path = "contract") {
    if (!is.list(value)) {
        return(invisible(TRUE))
    }
    forbidden <- c(
        "private_path",
        "source_path",
        "source_values",
        "source_rows",
        "source_identifiers",
        "source_annotations",
        "source_sequences",
        "source_distribution",
        "private_distribution",
        "observed_distribution",
        "empirical_quantiles",
        "private_records"
    )
    hits <- intersect(names(value), forbidden)
    if (length(hits)) {
        omicsWorkloadAbort(
            sprintf(
                "%s contains forbidden private-data fields: %s",
                path,
                paste(hits, collapse = ", ")
            ),
            class = "omics_workload_privacy_error"
        )
    }
    for (name in names(value)) {
        omicsWorkloadForbiddenPrivateFields(
            value[[name]],
            paste(path, name, sep = ".")
        )
    }
    invisible(TRUE)
}

omicsWorkloadValidateContract <- function(contract) {
    root_fields <- c(
        "schema",
        "schema_version",
        "workload_id",
        "capability",
        "adapter",
        "rng",
        "dimensions",
        "assay_mix",
        "distributions",
        "missingness",
        "censoring",
        "duplicate_policy",
        "sample_perturbations",
        "adapter_parameters",
        "privacy",
        "execution",
        "expected_digests"
    )
    omicsWorkloadRequireNames(contract, root_fields, "workload contract")
    if (!identical(contract$schema, .OMICS_WORKLOAD_SCHEMA)) {
        omicsWorkloadAbort(sprintf(
            "unsupported workload schema '%s'",
            as.character(contract$schema)
        ))
    }
    if (!identical(contract$schema_version, .OMICS_WORKLOAD_SCHEMA_VERSION)) {
        omicsWorkloadAbort(sprintf(
            "unsupported workload schema version '%s'",
            as.character(contract$schema_version)
        ))
    }
    omicsWorkloadRequireText(contract$workload_id, "workload_id")

    capability_fields <- c(
        "capability_id",
        "omic_type",
        "input_format",
        "data_level",
        "acquisition_mode"
    )
    omicsWorkloadRequireNames(
        contract$capability,
        capability_fields,
        "capability"
    )
    invisible(lapply(names(contract$capability), \(field) {
        omicsWorkloadRequireText(
            contract$capability[[field]],
            paste0("capability.", field)
        )
    }))

    omicsWorkloadRequireNames(
        contract$adapter,
        c("adapter_id", "adapter_version", "source_sha256"),
        "adapter"
    )
    omicsWorkloadRequireText(contract$adapter$adapter_id, "adapter.adapter_id")
    omicsWorkloadRequireText(
        contract$adapter$adapter_version,
        "adapter.adapter_version"
    )
    omicsWorkloadRequireDigest(
        contract$adapter$source_sha256,
        "adapter.source_sha256"
    )

    omicsWorkloadRequireNames(contract$rng, c("kind", "seed"), "rng")
    rng_kind <- unlist(contract$rng$kind, use.names = FALSE)
    if (!is.character(rng_kind) || length(rng_kind) != 3L ||
        anyNA(rng_kind) || any(!nzchar(rng_kind))) {
        omicsWorkloadAbort("rng.kind must contain the three RNGkind values")
    }
    omicsWorkloadRequireCount(contract$rng$seed, "rng.seed", minimum = 0L)

    omicsWorkloadRequireNames(
        contract$dimensions,
        c("feature_count", "sample_count", "assay_count"),
        "dimensions"
    )
    invisible(lapply(names(contract$dimensions), \(field) {
        omicsWorkloadRequireCount(
            contract$dimensions[[field]],
            paste0("dimensions.", field)
        )
    }))
    for (field in c(
        "assay_mix",
        "distributions",
        "missingness",
        "censoring",
        "duplicate_policy",
        "sample_perturbations",
        "adapter_parameters"
    )) {
        if (!is.list(contract[[field]])) {
            omicsWorkloadAbort(sprintf("%s must be a list", field))
        }
    }

    omicsWorkloadRequireNames(
        contract$privacy,
        c("classification", "source", "scale_metadata"),
        "privacy"
    )
    privacy_classes <- c("public_synthetic", "private_scale_only")
    if (!contract$privacy$classification %in% privacy_classes) {
        omicsWorkloadAbort(sprintf(
            "privacy.classification must be one of: %s",
            paste(privacy_classes, collapse = ", ")
        ))
    }
    omicsWorkloadRequireText(contract$privacy$source, "privacy.source")
    if (identical(contract$privacy$classification, "private_scale_only")) {
        omicsWorkloadRequireNames(
            contract$privacy$scale_metadata,
            c(
                "row_count",
                "column_count",
                "byte_size",
                "salted_source_fingerprint"
            ),
            "privacy.scale_metadata"
        )
        for (field in c("row_count", "column_count", "byte_size")) {
            omicsWorkloadRequireCount(
                contract$privacy$scale_metadata[[field]],
                paste0("privacy.scale_metadata.", field)
            )
        }
        omicsWorkloadRequireDigest(
            contract$privacy$scale_metadata$salted_source_fingerprint,
            "privacy.scale_metadata.salted_source_fingerprint"
        )
    } else if (!is.null(contract$privacy$scale_metadata)) {
        omicsWorkloadAbort(
            "privacy.scale_metadata must be null for public_synthetic workloads"
        )
    }

    execution_fields <- c(
        "repetitions",
        "sampling_interval_ms",
        "retention_window_ms",
        "timeout_ms",
        "maximum_retained_samples",
        "threads",
        "locale",
        "timezone",
        "cache_state"
    )
    omicsWorkloadRequireNames(contract$execution, execution_fields, "execution")
    for (field in c(
        "repetitions",
        "sampling_interval_ms",
        "retention_window_ms",
        "timeout_ms",
        "maximum_retained_samples",
        "threads"
    )) {
        omicsWorkloadRequireCount(
            contract$execution[[field]],
            paste0("execution.", field)
        )
    }
    invisible(lapply(c("locale", "timezone", "cache_state"), \(field) {
        omicsWorkloadRequireText(
            contract$execution[[field]],
            paste0("execution.", field)
        )
    }))

    omicsWorkloadRequireNames(
        contract$expected_digests,
        c("payload_sha256", "truth_sha256"),
        "expected_digests"
    )
    omicsWorkloadRequireDigest(
        contract$expected_digests$payload_sha256,
        "expected_digests.payload_sha256"
    )
    omicsWorkloadRequireDigest(
        contract$expected_digests$truth_sha256,
        "expected_digests.truth_sha256"
    )
    omicsWorkloadForbiddenPrivateFields(contract)
    invisible(contract)
}

omicsWorkloadCanonicalize <- function(value) {
    if (!is.list(value)) {
        return(value)
    }
    if (!is.null(names(value))) {
        value <- value[order(names(value), method = "radix")]
    }
    lapply(value, omicsWorkloadCanonicalize)
}

omicsWorkloadDigest <- function(value) {
    canonical <- omicsWorkloadCanonicalize(value)
    payload <- jsonlite::toJSON(
        canonical,
        auto_unbox = TRUE,
        null = "null",
        na = "null",
        digits = 17
    )
    digest::digest(payload, algo = "sha256", serialize = FALSE)
}

omicsWorkloadFileDigest <- function(path) {
    if (!file.exists(path) || dir.exists(path)) {
        omicsWorkloadAbort(sprintf("evidence file does not exist: %s", path))
    }
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

omicsWorkloadReadContract <- function(path) {
    if (!file.exists(path)) {
        omicsWorkloadAbort(sprintf("workload contract not found: %s", path))
    }
    contract <- jsonlite::read_json(path, simplifyVector = FALSE)
    omicsWorkloadValidateContract(contract)
    contract
}

omicsWorkloadLoadAdapter <- function(path, contract) {
    if (!file.exists(path)) {
        omicsWorkloadAbort(sprintf("workload adapter not found: %s", path))
    }
    observed_digest <- omicsWorkloadFileDigest(path)
    if (!identical(observed_digest, contract$adapter$source_sha256)) {
        omicsWorkloadAbort(
            "workload adapter source digest does not match the contract",
            class = "omics_workload_binding_error"
        )
    }
    adapter_environment <- new.env(parent = baseenv())
    sys.source(path, envir = adapter_environment)
    if (!exists("omicsWorkloadAdapter", envir = adapter_environment, inherits = FALSE)) {
        omicsWorkloadAbort("adapter must define omicsWorkloadAdapter()")
    }
    adapter <- adapter_environment$omicsWorkloadAdapter()
    omicsWorkloadRequireNames(
        adapter,
        c("adapter_id", "adapter_version", "supported_omics", "prepare", "run"),
        "adapter result"
    )
    omicsWorkloadRequireText(adapter$adapter_id, "adapter result.adapter_id")
    omicsWorkloadRequireText(
        adapter$adapter_version,
        "adapter result.adapter_version"
    )
    supported_omics <- unlist(adapter$supported_omics, use.names = FALSE)
    if (!is.character(supported_omics) || !length(supported_omics)) {
        omicsWorkloadAbort("adapter result.supported_omics must contain omic names")
    }
    if (!is.function(adapter$prepare) || !is.function(adapter$run)) {
        omicsWorkloadAbort("adapter prepare and run entries must be functions")
    }
    matches <- identical(adapter$adapter_id, contract$adapter$adapter_id) &&
        identical(adapter$adapter_version, contract$adapter$adapter_version) &&
        contract$capability$omic_type %in% supported_omics
    if (!matches) {
        omicsWorkloadAbort(
            "adapter identity or supported omic does not match the contract",
            class = "omics_workload_binding_error"
        )
    }
    adapter
}

omicsWorkloadEvidenceBinding <- function(
    contract,
    adapter_path,
    payload_path,
    truth_path,
    code_revision,
    environment
) {
    omicsWorkloadValidateContract(contract)
    binding <- list(
        schema = .OMICS_WORKLOAD_EVIDENCE_SCHEMA,
        schema_version = .OMICS_WORKLOAD_EVIDENCE_VERSION,
        workload_id = contract$workload_id,
        capability = contract$capability,
        contract_sha256 = omicsWorkloadDigest(contract),
        adapter = c(
            contract$adapter,
            list(observed_source_sha256 = omicsWorkloadFileDigest(adapter_path))
        ),
        payload_sha256 = omicsWorkloadFileDigest(payload_path),
        truth_sha256 = omicsWorkloadFileDigest(truth_path),
        code_revision = code_revision,
        environment = environment,
        privacy_classification = contract$privacy$classification
    )
    expected <- contract$expected_digests
    if (!identical(binding$payload_sha256, expected$payload_sha256) ||
        !identical(binding$truth_sha256, expected$truth_sha256)) {
        omicsWorkloadAbort(
            "generated payload or truth digest does not match the contract",
            class = "omics_workload_binding_error"
        )
    }
    binding$binding_sha256 <- omicsWorkloadDigest(binding)
    binding
}
