# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

.LIPID_CLOSEOUT_SCHEMA <- "multischolar.lipidomics_artifact_closeout"
.LIPID_CLOSEOUT_VERSION <- "1.0.0"
.LIPID_CLOSEOUT_REQUIRED_PAIRS <- 5L

artifactLipidomicsPromotionGates <- function() {
    list(
        minimum_peak_rss_reduction_fraction = 0.30,
        minimum_retained_rss_reduction_fraction = 0.40,
        maximum_runtime_ratio = 1.25,
        maximum_committed_disk_ratio = 1.35,
        maximum_peak_disk_ratio = 2.20,
        maximum_final_file_count_ratio = 4.00,
        maximum_bounded_query_p95_seconds = 0.50,
        minimum_representative_repetitions =
            .LIPID_CLOSEOUT_REQUIRED_PAIRS
    )
}

artifactLipidomicsCloseoutEvidence <- function() {
    list(
        workload_id = paste0(
            "lipidomics.custom.mixed.public-representative-",
            "16000x12.v1"
        ),
        generator_version = "1.0.0",
        assay_order = c("LCMS_Pos", "LCMS_Neg", "GCMS"),
        dimensions = list(features = 16000L, samples = 12L, assays = 3L),
        paired_repetitions = 5L,
        runtime_paired_repetitions = 1L,
        required_paired_repetitions = .LIPID_CLOSEOUT_REQUIRED_PAIRS,
        memory = list(
            settled_rss_ratio = 1,
            elapsed_seconds = 0.302
        ),
        artifact = list(
            settled_rss_ratio = 1.0212,
            elapsed_seconds = 5.66,
            source_payloads_retained = TRUE,
            eviction_performed = FALSE
        ),
        rss_campaign = list(
            profiles = c(
                lcms_positive = 1.0325,
                lcms_negative = 1.0331,
                gcms = 1.0349,
                mixed = 1.0212,
                private_calibrated = 1.0221
            ),
            zero_reduction_runs = 23L,
            total_artifact_runs = 25L,
            best_reduction_fractions = c(0.0098, 0.0104)
        ),
        eviction_gate = list(
            ticket_id = "OMICS-ART-052",
            status = "blocked",
            reason = "representative_rss_benefit_not_demonstrated"
        ),
        private_evidence = list(
            pooled_with_public = FALSE,
            source_path_retained = FALSE,
            unsalted_fingerprint_retained = FALSE,
            headers_or_values_retained = FALSE
        ),
        comparable = TRUE,
        claim_boundary = paste(
            "generated workload authorizes operational performance, scale,",
            "query, portability, and lifecycle evidence only"
        )
    )
}

artifactLipidomicsSupportedCloseout <- function(descriptor) {
    evidence <- artifactLipidomicsCloseoutEvidence()
    gates <- artifactLipidomicsPromotionGates()
    rss_reduction <- 1 - (
        evidence$artifact$settled_rss_ratio /
            evidence$memory$settled_rss_ratio
    )
    runtime_ratio <- evidence$artifact$elapsed_seconds /
        evidence$memory$elapsed_seconds
    checks <- c(
        enough_pairs = evidence$paired_repetitions >=
            gates$minimum_representative_repetitions,
        enough_runtime_pairs = evidence$runtime_paired_repetitions >=
            gates$minimum_representative_repetitions,
        eviction_passed = identical(evidence$eviction_gate$status, "passed"),
        peak_rss = rss_reduction >=
            gates$minimum_peak_rss_reduction_fraction,
        runtime = runtime_ratio <= gates$maximum_runtime_ratio,
        finite = all(is.finite(c(
            evidence$memory$settled_rss_ratio,
            evidence$memory$elapsed_seconds,
            evidence$artifact$settled_rss_ratio,
            evidence$artifact$elapsed_seconds
        ))),
        comparable = isTRUE(evidence$comparable)
    )
    list(
        schema = .LIPID_CLOSEOUT_SCHEMA,
        schema_version = .LIPID_CLOSEOUT_VERSION,
        capability_id = descriptor$descriptor_id,
        capability_version = descriptor$descriptor_version,
        descriptor_digest = descriptor$descriptor_digest,
        support_status = "scientifically_supported",
        scientific_evidence = list(
            fixture_ids = descriptor$fixtures$fixture_ids,
            oracle_id = descriptor$scientific_oracle$oracle_id,
            formats = "lipidsearch"
        ),
        performance_evidence = c(evidence, list(
            gates = gates,
            checks = checks,
            measured = list(
                peak_rss_reduction_fraction = rss_reduction,
                runtime_ratio = runtime_ratio
            )
        )),
        certification = descriptor$certification,
        promotion_status = "withheld",
        reason_code = "eviction_and_performance_gates_failed",
        effective_default_backend = "memory",
        maximum_forced_rollout = "dual_write",
        existing_projects_require_explicit_migration = TRUE,
        rollback = descriptor$rollback
    )
}

artifactLipidomicsUnsupportedCloseout <- function(capability) {
    list(
        schema = .LIPID_CLOSEOUT_SCHEMA,
        schema_version = .LIPID_CLOSEOUT_VERSION,
        capability_id = capability$capability_id,
        capability_version = capability$capability_version,
        descriptor_digest = NULL,
        support_status = "reader_characterized",
        scientific_evidence = list(
            fixture_ids = list(),
            oracle_id = NULL,
            formats = capability$identity$input_format
        ),
        performance_evidence = list(
            paired_repetitions = 0L,
            required_paired_repetitions = .LIPID_CLOSEOUT_REQUIRED_PAIRS,
            gates = artifactLipidomicsPromotionGates(),
            checks = c(enough_pairs = FALSE)
        ),
        certification = list(status = "uncertified", auto_eligible = FALSE),
        promotion_status = "withheld",
        reason_code = "complete_scientific_workflow_unverified",
        effective_default_backend = "memory",
        maximum_forced_rollout = "none",
        existing_projects_require_explicit_migration = TRUE,
        rollback = list(
            owner_id = capability$capability_id,
            strategy_id = "retain_memory_backend",
            target_backend = "memory"
        )
    )
}

validateArtifactLipidomicsCloseout <- function(decision) {
    required <- c(
        "schema", "schema_version", "capability_id", "capability_version",
        "descriptor_digest", "support_status", "scientific_evidence",
        "performance_evidence", "certification", "promotion_status",
        "reason_code", "effective_default_backend", "maximum_forced_rollout",
        "existing_projects_require_explicit_migration", "rollback"
    )
    valid <- is.list(decision) && identical(names(decision), required) &&
        identical(decision$schema, .LIPID_CLOSEOUT_SCHEMA) &&
        identical(decision$schema_version, .LIPID_CLOSEOUT_VERSION) &&
        workflowCapabilityScalarString(decision$capability_id) &&
        workflowCapabilityScalarString(decision$capability_version) &&
        decision$support_status %in% c(
            "scientifically_supported", "reader_characterized"
        ) && identical(decision$promotion_status, "withheld") &&
        identical(decision$effective_default_backend, "memory") &&
        identical(decision$existing_projects_require_explicit_migration, TRUE)
    if (!isTRUE(valid)) {
        workflowSessionAbort(
            "lipidomics closeout decision is malformed or promotes unsafely",
            "multischolar_invalid_lipid_closeout"
        )
    }
    if (!is.null(decision$descriptor_digest)) {
        artifactRefValidateDigest(
            decision$descriptor_digest,
            "lipidomics closeout descriptor digest"
        )
    }
    checks <- decision$performance_evidence$checks
    if (identical(decision$support_status, "scientifically_supported") &&
        (is.null(checks) || all(checks))) {
        workflowSessionAbort(
            "lipidomics closeout evidence unexpectedly satisfies promotion",
            "multischolar_incomplete_lipid_closeout_evidence"
        )
    }
    decision
}

artifactLipidomicsCloseoutDecisions <- function() {
    descriptor <- artifactLipidomicsWorkflowDescriptor()
    capabilities <- workflowCapabilityCatalogue()
    custom <- artifactLipidomicsSupportedCloseout(descriptor)
    others <- Filter(function(capability) {
        identical(capability$identity$omic_type, "lipidomics") &&
            !identical(capability$capability_id, descriptor$descriptor_id)
    }, capabilities)
    decisions <- c(
        stats::setNames(list(custom), custom$capability_id),
        lapply(others, artifactLipidomicsUnsupportedCloseout)
    )
    decisions <- lapply(decisions, validateArtifactLipidomicsCloseout)
    names(decisions) <- vapply(
        decisions,
        `[[`,
        character(1),
        "capability_id"
    )
    decisions
}
