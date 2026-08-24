# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

.PROT_NONDIA_CLOSEOUT_SCHEMA <-
    "multischolar.proteomics_nondia_artifact_closeout"
.PROT_NONDIA_CLOSEOUT_VERSION <- "1.0.0"
.PROT_NONDIA_CLOSEOUT_REQUIRED_PAIRS <- 5L

#' Define unchanged non-DIA promotion performance gates
#'
#' @return A named list of immutable paired-performance gates.
#' @noRd
artifactProteomicsNonDiaPromotionGates <- function() {
    list(
        minimum_peak_rss_reduction_fraction = 0.30,
        minimum_retained_rss_reduction_fraction = 0.40,
        maximum_runtime_ratio = 1.25,
        maximum_committed_disk_ratio = 1.35,
        maximum_peak_disk_ratio = 2.20,
        maximum_final_file_count_ratio = 4.00,
        maximum_bounded_query_p95_ratio = 2.00,
        maximum_bounded_query_p95_seconds = 0.50,
        minimum_representative_repetitions =
            .PROT_NONDIA_CLOSEOUT_REQUIRED_PAIRS
    )
}

#' Resolve the maximum rollout declared by one descriptor
#'
#' @param descriptor Validated workflow descriptor.
#'
#' @return The highest declared rollout stage.
#' @noRd
artifactProteomicsNonDiaMaximumRollout <- function(descriptor) {
    descriptor <- validateArtifactWorkflowDescriptor(descriptor)
    rollouts <- vapply(
        descriptor$stages,
        `[[`,
        character(1),
        "maximum_rollout"
    )
    .WORKFLOW_ARTIFACT_ROLLOUTS[[max(match(
        rollouts,
        .WORKFLOW_ARTIFACT_ROLLOUTS
    ))]]
}

#' Build one supported non-DIA closeout decision
#'
#' @param descriptor Validated non-DIA workflow descriptor.
#' @param spec Tuple-specific evidence specification.
#'
#' @return A data-only fail-closed closeout decision.
#' @noRd
artifactProteomicsNonDiaSupportedCloseout <- function(descriptor, spec) {
    list(
        schema = .PROT_NONDIA_CLOSEOUT_SCHEMA,
        schema_version = .PROT_NONDIA_CLOSEOUT_VERSION,
        capability_id = descriptor$descriptor_id,
        capability_version = descriptor$descriptor_version,
        descriptor_digest = descriptor$descriptor_digest,
        support_status = "scientifically_supported",
        scientific_evidence = list(
            fixture_ids = descriptor$fixtures$fixture_ids,
            oracle_id = descriptor$scientific_oracle$oracle_id,
            e2e_lane = spec$workflow_slug
        ),
        performance_evidence = list(
            memory_baseline = list(
                status = "passed_historical_memory_only",
                evidence_id = spec$resource_id,
                repetitions = 5L
            ),
            paired_artifact = list(
                status = "absent",
                valid_pairs = 0L,
                required_pairs = .PROT_NONDIA_CLOSEOUT_REQUIRED_PAIRS
            ),
            gates = artifactProteomicsNonDiaPromotionGates()
        ),
        certification = list(
            status = descriptor$certification$status,
            auto_eligible = descriptor$certification$auto_eligible
        ),
        maximum_rollout = artifactProteomicsNonDiaMaximumRollout(descriptor),
        promotion_status = "withheld",
        reason_code = "paired_artifact_performance_evidence_absent",
        effective_default_backend = "memory",
        forced_artifact_canary = TRUE,
        rollback = descriptor$rollback
    )
}

#' Build one unsupported non-DIA closeout decision
#'
#' @param capability Unverified workflow capability.
#'
#' @return A data-only unsupported closeout decision.
#' @noRd
artifactProteomicsNonDiaUnsupportedCloseout <- function(capability) {
    list(
        schema = .PROT_NONDIA_CLOSEOUT_SCHEMA,
        schema_version = .PROT_NONDIA_CLOSEOUT_VERSION,
        capability_id = capability$capability_id,
        capability_version = capability$capability_version,
        descriptor_digest = NULL,
        support_status = "advertised_unverified",
        scientific_evidence = list(
            fixture_ids = list(),
            oracle_id = NULL,
            e2e_lane = NULL
        ),
        performance_evidence = list(
            memory_baseline = list(
                status = "absent",
                evidence_id = NULL,
                repetitions = 0L
            ),
            paired_artifact = list(
                status = "absent",
                valid_pairs = 0L,
                required_pairs = .PROT_NONDIA_CLOSEOUT_REQUIRED_PAIRS
            ),
            gates = artifactProteomicsNonDiaPromotionGates()
        ),
        certification = list(status = "uncertified", auto_eligible = FALSE),
        maximum_rollout = "none",
        promotion_status = "withheld",
        reason_code = "scientific_support_evidence_absent",
        effective_default_backend = "memory",
        forced_artifact_canary = FALSE,
        rollback = list(
            owner_id = capability$capability_id,
            strategy_id = "retain_memory_backend",
            target_backend = "memory"
        )
    )
}

#' Validate one non-DIA closeout decision
#'
#' @param decision Candidate closeout decision.
#'
#' @return The validated decision.
#' @noRd
validateArtifactProteomicsNonDiaCloseout <- function(decision) {
    required <- c(
        "schema", "schema_version", "capability_id", "capability_version",
        "descriptor_digest", "support_status", "scientific_evidence",
        "performance_evidence", "certification", "maximum_rollout",
        "promotion_status", "reason_code", "effective_default_backend",
        "forced_artifact_canary", "rollback"
    )
    valid <- is.list(decision) && identical(names(decision), required) &&
        identical(decision$schema, .PROT_NONDIA_CLOSEOUT_SCHEMA) &&
        identical(decision$schema_version, .PROT_NONDIA_CLOSEOUT_VERSION) &&
        workflowCapabilityScalarString(decision$capability_id) &&
        workflowCapabilityScalarString(decision$capability_version) &&
        decision$support_status %in% c(
            "scientifically_supported",
            "advertised_unverified"
        ) && identical(decision$promotion_status, "withheld") &&
        identical(decision$effective_default_backend, "memory") &&
        is.logical(decision$forced_artifact_canary) &&
        length(decision$forced_artifact_canary) == 1L &&
        is.list(decision$scientific_evidence) &&
        is.list(decision$performance_evidence) &&
        is.list(decision$certification) && is.list(decision$rollback)
    if (!isTRUE(valid)) {
        artifactDescriptorAbort(
            "non-DIA proteomics closeout decision is malformed",
            "multischolar_invalid_prot_nondia_closeout"
        )
    }
    paired <- decision$performance_evidence$paired_artifact
    withheld <- identical(paired$status, "absent") &&
        identical(as.integer(paired$valid_pairs), 0L) &&
        identical(
            as.integer(paired$required_pairs),
            .PROT_NONDIA_CLOSEOUT_REQUIRED_PAIRS
        ) && !isTRUE(decision$certification$auto_eligible)
    if (!isTRUE(withheld)) {
        artifactDescriptorAbort(
            "non-DIA closeout cannot promote without paired artifact evidence",
            "multischolar_incomplete_prot_nondia_closeout_evidence"
        )
    }
    if (!is.null(decision$descriptor_digest)) {
        artifactRefValidateDigest(
            decision$descriptor_digest,
            "closeout descriptor digest"
        )
    }
    decision
}

#' Build all exact non-DIA proteomics closeout decisions
#'
#' @return Five independently validated decisions keyed by capability id.
#' @noRd
artifactProteomicsNonDiaCloseoutDecisions <- function() {
    descriptors <- artifactProteomicsNonDiaWorkflowDescriptors()
    specs <- artifactProteomicsNonDiaDescriptorSpecs()
    supported <- Map(
        artifactProteomicsNonDiaSupportedCloseout,
        descriptors,
        specs[names(descriptors)]
    )
    capabilities <- workflowCapabilityCatalogue()
    unsupported <- Filter(
        \(capability) identical(capability$identity$input_format, "spectronaut"),
        capabilities
    )
    unsupported <- lapply(
        unsupported,
        artifactProteomicsNonDiaUnsupportedCloseout
    )
    decisions <- c(supported, unsupported)
    decisions <- lapply(decisions, validateArtifactProteomicsNonDiaCloseout)
    stats::setNames(
        decisions,
        vapply(decisions, `[[`, character(1), "capability_id")
    )
}
