# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

protDiaDaValidateParameters <- function(parameters) {
    required <- c(
        "formula_string", "da_q_val_thresh", "treat_lfc_cutoff",
        "eBayes_trend", "eBayes_robust", "qvalue_column",
        "raw_pvalue_column", "fdr_fallback"
    )
    numeric_valid <- function(value) {
        length(value) == 1L && is.numeric(value) && !is.na(value) &&
            is.finite(value)
    }
    valid <- is.list(parameters) && identical(names(parameters), required) &&
        workflowCapabilityScalarString(parameters$formula_string) &&
        numeric_valid(parameters$da_q_val_thresh) &&
        parameters$da_q_val_thresh >= 0 && parameters$da_q_val_thresh <= 1 &&
        numeric_valid(parameters$treat_lfc_cutoff) &&
        parameters$treat_lfc_cutoff >= 0 &&
        all(vapply(
            parameters[c("eBayes_trend", "eBayes_robust")],
            function(value) is.logical(value) && length(value) == 1L &&
                !is.na(value),
            logical(1)
        )) && identical(parameters$qvalue_column, "fdr_qvalue") &&
        identical(parameters$raw_pvalue_column, "raw_pvalue") &&
        identical(parameters$fdr_fallback, "BH")
    if (!isTRUE(valid)) {
        protDiaDaArtifactAbort(
            "DIA-NN DA parameters are malformed or unsupported",
            "multischolar_invalid_prot_dia_da_manifest"
        )
    }
    parameters$da_q_val_thresh <- as.double(parameters$da_q_val_thresh)
    parameters$treat_lfc_cutoff <- as.double(parameters$treat_lfc_cutoff)
    parameters
}

protDiaDaValidateSoftware <- function(software) {
    required <- c("multischolar", "limma", "qvalue", "r")
    valid <- is.list(software) && identical(names(software), required) &&
        all(vapply(software, workflowCapabilityScalarString, logical(1)))
    if (!isTRUE(valid)) {
        protDiaDaArtifactAbort(
            "DIA-NN DA software provenance is malformed",
            "multischolar_invalid_prot_dia_da_manifest"
        )
    }
    software
}

protDiaDaValidateSummary <- function(summary) {
    required <- c(
        "rows", "significant", "up", "down", "q_value_threshold",
        "fold_change_cutoff"
    )
    counts <- unlist(summary[c("rows", "significant", "up", "down")])
    thresholds <- unlist(summary[c("q_value_threshold", "fold_change_cutoff")])
    valid <- is.list(summary) && identical(names(summary), required) &&
        is.numeric(counts) && length(counts) == 4L && !anyNA(counts) &&
        all(is.finite(counts)) && all(counts >= 0) &&
        all(counts == floor(counts)) && all(counts <= counts[[1L]]) &&
        is.numeric(thresholds) && length(thresholds) == 2L &&
        !anyNA(thresholds) && all(is.finite(thresholds)) &&
        thresholds[[1L]] >= 0 && thresholds[[1L]] <= 1 &&
        thresholds[[2L]] >= 0
    if (!isTRUE(valid)) {
        protDiaDaArtifactAbort(
            "DIA-NN DA summary is malformed",
            "multischolar_invalid_prot_dia_da_manifest"
        )
    }
    summary$rows <- as.integer(summary$rows)
    summary$significant <- as.integer(summary$significant)
    summary$up <- as.integer(summary$up)
    summary$down <- as.integer(summary$down)
    summary$q_value_threshold <- as.double(summary$q_value_threshold)
    summary$fold_change_cutoff <- as.double(summary$fold_change_cutoff)
    summary
}

protDiaDaValidateRunEntry <- function(entry) {
    required <- c(
        "contrast_id", "contrast_name", "full_format", "friendly_name",
        "manifest_relative_path", "manifest_digest"
    )
    valid <- is.list(entry) && identical(names(entry), required) &&
        all(vapply(entry, workflowCapabilityScalarString, logical(1)))
    if (!isTRUE(valid)) {
        protDiaDaArtifactAbort(
            "DIA-NN DA run contrast entry is malformed",
            "multischolar_invalid_prot_dia_da_manifest"
        )
    }
    artifactValidatePathComponent(entry$contrast_id, "da_contrast_id")
    artifactNormalizeRelativePath(entry$manifest_relative_path)
    artifactRefValidateDigest(entry$manifest_digest, "contrast manifest digest")
    entry
}

protDiaDaValidateRunManifest <- function(manifest, context = NULL) {
    required <- c(
        "schema", "schema_version", "project_id", "workflow_id",
        "descriptor_contract", "run_id", "source", "parameters",
        "parameters_digest", "software", "contrasts", "created_at",
        "manifest_digest"
    )
    valid <- is.list(manifest) && identical(names(manifest), required) &&
        identical(manifest$schema, .PROT_DIA_DA_RUN_SCHEMA) &&
        identical(
            workflowStateVersionValue(manifest$schema_version),
            .PROT_DIA_DA_RUN_VERSION
        ) && all(vapply(
            manifest[c("project_id", "workflow_id", "run_id")],
            workflowCapabilityScalarString,
            logical(1)
        )) && is.list(manifest$source) && is.list(manifest$contrasts) &&
        length(manifest$contrasts) > 0L && artifactRefValidUtc(manifest$created_at)
    if (!isTRUE(valid)) {
        protDiaDaArtifactAbort(
            "DIA-NN DA run manifest is malformed or unsupported",
            "multischolar_invalid_prot_dia_da_manifest"
        )
    }
    artifactWorkflowStateAssertSafeMetadata(manifest, "DIA-NN DA run manifest")
    if (!identical(manifest$manifest_digest, protDiaDaJsonDigest(manifest))) {
        protDiaDaArtifactAbort(
            "DIA-NN DA run manifest fingerprint differs",
            "multischolar_prot_dia_da_manifest_digest_mismatch"
        )
    }
    parameters <- protDiaDaValidateParameters(manifest$parameters)
    software <- protDiaDaValidateSoftware(manifest$software)
    entries <- lapply(manifest$contrasts, protDiaDaValidateRunEntry)
    entry_ids <- vapply(entries, `[[`, character(1), "contrast_id")
    entry_names <- vapply(entries, `[[`, character(1), "contrast_name")
    if (anyDuplicated(entry_ids) > 0L || anyDuplicated(entry_names) > 0L ||
        !identical(
            manifest$parameters_digest,
            artifactSemanticDigest(parameters)
        )) {
        protDiaDaArtifactAbort(
            "DIA-NN DA run manifest bindings differ",
            "multischolar_prot_dia_da_manifest_digest_mismatch"
        )
    }
    contrasts <- protDiaDaRunContrasts(data.frame(
        contrasts = entry_names,
        full_format = vapply(entries, `[[`, character(1), "full_format"),
        friendly_names = vapply(entries, `[[`, character(1), "friendly_name")
    ))
    expected_run_id <- protDiaDaRunId(
        manifest$source, contrasts, parameters, software
    )
    if (!identical(manifest$run_id, expected_run_id)) {
        protDiaDaArtifactAbort(
            "DIA-NN DA run ID differs from its scientific inputs",
            "multischolar_prot_dia_da_manifest_digest_mismatch"
        )
    }
    if (!is.null(context)) {
        protDiaDaValidateRunContext(manifest, entries, context)
    }
    manifest$schema_version <- .PROT_DIA_DA_RUN_VERSION
    manifest$parameters <- parameters
    manifest$software <- software
    manifest$contrasts <- entries
    manifest
}

protDiaDaValidateRunContext <- function(manifest, entries, context) {
    identity <- context$getIdentity()
    descriptor <- findArtifactWorkflowDescriptor(
        identity, artifactWorkflowDescriptorCatalogue()
    )
    expected_contract <- list(
        descriptor_id = descriptor$descriptor_id,
        descriptor_version = descriptor$descriptor_version,
        descriptor_digest = descriptor$descriptor_digest
    )
    if (!identical(manifest$project_id, identity$project_id) ||
        !identical(manifest$workflow_id, identity$workflow_id) ||
        !identical(manifest$descriptor_contract, expected_contract)) {
        protDiaDaArtifactAbort(
            "DIA-NN DA run belongs to a different workflow contract",
            "multischolar_prot_dia_da_identity_mismatch"
        )
    }
    store <- protDiaDaStore(context)
    protDiaDaValidateRunSource(manifest$source, store)
    invisible(lapply(entries, function(entry) {
        path <- artifactResolveContainedPath(
            context$getProjectRoot(),
            entry$manifest_relative_path,
            must_exist = TRUE
        )
        contrast <- protDiaDaReadJson(
            path,
            function(value) protDiaDaValidateContrastManifest(value, store)
        )
        binding <- contrast$contrast
        names(binding) <- c(
            "contrast_id", "contrast_name", "full_format", "friendly_name"
        )
        valid <- identical(contrast$manifest_digest, entry$manifest_digest) &&
            identical(contrast$run_id, manifest$run_id) &&
            identical(contrast$project_id, manifest$project_id) &&
            identical(contrast$workflow_id, manifest$workflow_id) &&
            identical(
                contrast$source_generation_id,
                manifest$source$generation_id
            ) && identical(binding, entry[names(binding)])
        if (!isTRUE(valid)) {
            protDiaDaArtifactAbort(
                "DIA-NN DA run contrast binding differs",
                "multischolar_prot_dia_da_contrast_mismatch"
            )
        }
        invisible(entry)
    }))
    invisible(TRUE)
}
