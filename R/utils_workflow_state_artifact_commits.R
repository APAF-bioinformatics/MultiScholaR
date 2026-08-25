# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

artifactWorkflowStateRepresentation <- function(persistence_hint) {
    if (is.null(persistence_hint)) return("materialized")
    representation <- persistence_hint$representation
    if (!workflowStateScalarString(representation) ||
        !representation %in% c(
            "row_selection", "assay_selection", "state_reference"
        )) {
        artifactWorkflowStateAbort(
            "artifact persistence hint has an unsupported representation",
            "multischolar_invalid_artifact_state_persistence_hint"
        )
    }
    representation
}

artifactWorkflowStateDataMetrics <- function(store, data, state_object) {
    reference_names <- names(data$artifact_refs)
    reference_bytes <- vapply(data$artifact_refs, \(ref) {
        path <- artifactStoreResolveFile(
            store,
            ref$relative_path,
            must_exist = TRUE
        )
        unname(as.numeric(file.info(path)$size))
    }, numeric(1))
    new_reference_names <- data$new_reference_names
    if (is.null(new_reference_names)) new_reference_names <- reference_names
    new_reference_names <- intersect(new_reference_names, reference_names)
    list(
        hydrated_object_bytes = unname(as.numeric(utils::object.size(state_object))),
        generation_artifact_bytes = unname(sum(reference_bytes)),
        new_artifact_bytes = unname(sum(reference_bytes[new_reference_names])),
        generation_artifact_count = as.integer(length(reference_names)),
        new_artifact_count = as.integer(length(new_reference_names))
    )
}

artifactWorkflowStateInvokeRegistryFailure <- function(
    failure_injector,
    generation_id,
    state_name,
    parent_generation_id
) {
    artifactStoreInvokeFailure(
        failure_injector,
        "before_state_registry_commit",
        list(
            generation_id = generation_id,
            state_name = state_name,
            parent_generation_id = parent_generation_id
        )
    )
}

artifactWorkflowStateRegisterData <- function(session, identity, data, timestamp) {
    if (isTRUE(data$reused)) return(invisible(FALSE))
    new_reference_names <- data$new_reference_names
    if (is.null(new_reference_names)) {
        new_reference_names <- names(data$artifact_refs)
    }
    for (index in seq_along(new_reference_names)) {
        reference_name <- new_reference_names[[index]]
        projectRegistryWrite(
            session,
            "artifact",
            artifactWorkflowStateArtifactRecord(
                identity,
                data$artifact_refs[[reference_name]],
                index - 1L
            )
        )
    }
    dependencies <- data$dependencies
    if (is.null(dependencies)) dependencies <- list()
    for (dependency in dependencies) {
        projectRegistryWrite(
            session,
            "dependency",
            c(
                list(workflow_id = identity$workflow_id),
                dependency,
                list(recorded_at = timestamp)
            )
        )
    }
    invisible(TRUE)
}
