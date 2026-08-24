# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

#' Test artifact session eligibility by exact proteomics descriptor
#' @param workflow_data Mutable proteomics workflow state.
#' @param direction Session export or restore direction.
#' @return A scalar logical.
#' @noRd
protArtifactSessionEligible <- function(workflow_data, direction = "export") {
    protDiaSessionArtifactEligible(workflow_data, direction) ||
        protNonDiaSessionArtifactEligible(workflow_data, direction)
}

#' Dispatch portable proteomics session publication
#' @param workflow_data Mutable proteomics workflow state.
#' @param norm_data Mutable normalization module state.
#' @param session_data Collected legacy-compatible session data.
#' @param export_artifacts Legacy export paths.
#' @param source_dir Project source directory.
#' @param failure_injector Optional session failure injector used by tests.
#' @return Published session manifest paths and data.
#' @noRd
saveProtArtifactSessionManifest <- function(
    workflow_data,
    norm_data,
    session_data,
    export_artifacts,
    source_dir,
    failure_injector = NULL
) {
    if (protDiaSessionArtifactEligible(workflow_data, "export")) {
        return(saveProtDiaSessionManifest(
            workflow_data,
            norm_data,
            session_data,
            export_artifacts,
            source_dir,
            failure_injector
        ))
    }
    saveProtNonDiaSessionManifest(
        workflow_data,
        norm_data,
        session_data,
        export_artifacts,
        source_dir,
        failure_injector
    )
}

#' Restore a portable proteomics session by descriptor ID
#' @param manifest_path Session manifest path.
#' @param experiment_paths Moved workflow project paths.
#' @param resource_policy Optional project registry resource policy.
#' @return A verified artifact session bundle.
#' @noRd
restoreProtArtifactSessionManifest <- function(
    manifest_path,
    experiment_paths,
    resource_policy = NULL
) {
    manifest <- readWorkflowSessionManifest(manifest_path)
    descriptor_id <- manifest$descriptor$descriptor_id
    if (identical(
        descriptor_id,
        artifactDiaWorkflowDescriptor()$descriptor_id
    )) {
        if (!identical(protDiaSessionMode("restore"), "enabled")) {
            workflowSessionAbort(
                "DIA-NN artifact session restore is disabled",
                "multischolar_prot_dia_session_restore_disabled"
            )
        }
        return(restoreProtDiaSessionManifest(
            manifest_path,
            experiment_paths,
            resource_policy
        ))
    }
    protNonDiaReadthroughDescriptor(descriptor_id)
    if (!identical(protNonDiaSessionMode("restore"), "enabled")) {
        workflowSessionAbort(
            "non-DIA artifact session restore is disabled",
            "multischolar_prot_nondia_session_restore_disabled"
        )
    }
    restoreProtNonDiaSessionManifest(
        manifest_path,
        experiment_paths,
        resource_policy
    )
}

#' Test whether the exact artifact session restore switch is enabled
#' @param manifest_path Candidate latest artifact manifest path.
#' @return A scalar logical; malformed manifests remain candidates for rejection.
#' @noRd
protArtifactSessionRestoreEnabled <- function(manifest_path) {
    if (!file.exists(manifest_path)) return(FALSE)
    manifest <- tryCatch(
        readWorkflowSessionManifest(manifest_path),
        error = \(error) error
    )
    if (inherits(manifest, "error")) {
        return(identical(protDiaSessionMode("restore"), "enabled"))
    }
    if (identical(
        manifest$descriptor$descriptor_id,
        artifactDiaWorkflowDescriptor()$descriptor_id
    )) {
        return(identical(protDiaSessionMode("restore"), "enabled"))
    }
    identical(protNonDiaSessionMode("restore"), "enabled")
}

#' Write a legacy compatibility session by exact descriptor
#' @param bundle Verified artifact session bundle.
#' @param path Destination RDS path.
#' @return Compatibility path and integrity fingerprint.
#' @noRd
writeProtArtifactCompatibilitySession <- function(bundle, path) {
    descriptor_id <- bundle$manifest$descriptor$descriptor_id
    if (identical(
        descriptor_id,
        artifactDiaWorkflowDescriptor()$descriptor_id
    )) {
        return(writeProtDiaCompatibilitySession(bundle, path))
    }
    protNonDiaReadthroughDescriptor(descriptor_id)
    writeProtNonDiaCompatibilitySession(bundle, path)
}
