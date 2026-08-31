# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

.WORKFLOW_CAPABILITY_KEY_FIELDS <- c(
    "omic_type",
    "workflow_id",
    "workflow_type",
    "workflow_slug",
    "input_format",
    "data_level",
    "acquisition_mode"
)
.WORKFLOW_IDENTITY_FIELDS <- c(
    "project_id",
    "omic_type",
    "omic_label",
    "workflow_id",
    "workflow_type",
    "workflow_slug",
    "input_format",
    "data_level",
    "acquisition_mode"
)
.WORKFLOW_BACKENDS <- c("memory", "artifact", "auto")
.WORKFLOW_ARTIFACT_ROLLOUTS <- c("dual_write", "read_through", "evict")
.WORKFLOW_PROJECT_STATES <- c(
    "new",
    "legacy_memory",
    "artifact_valid",
    "artifact_corrupt",
    "artifact_future_schema"
)

workflowCapabilityAbort <- function(message, class) {
    rlang::abort(message, class = c(class, "multischolar_workflow_error"))
}

workflowCapabilityScalarString <- function(value) {
    is.character(value) && length(value) == 1L && !is.na(value) && nzchar(value)
}

workflowCapabilityKey <- function(identity) {
    missing_fields <- setdiff(.WORKFLOW_CAPABILITY_KEY_FIELDS, names(identity))
    if (length(missing_fields) > 0L ||
        !all(vapply(
            identity[.WORKFLOW_CAPABILITY_KEY_FIELDS],
            workflowCapabilityScalarString,
            logical(1)
        ))) {
        workflowCapabilityAbort(
            "workflow identity does not contain a complete capability key",
            "multischolar_invalid_workflow_identity"
        )
    }
    paste(
        unlist(identity[.WORKFLOW_CAPABILITY_KEY_FIELDS], use.names = FALSE),
        collapse = "\x1f"
    )
}

newWorkflowCapability <- function(
    capability_id,
    identity,
    artifact_eligible = FALSE,
    auto_eligible = FALSE,
    maximum_artifact_rollout = NULL,
    explicit_maximum_artifact_rollout = maximum_artifact_rollout,
    capability_version = "1.0.0"
) {
    list(
        capability_id = capability_id,
        capability_version = capability_version,
        identity = identity,
        artifact_eligible = isTRUE(artifact_eligible),
        auto_eligible = isTRUE(auto_eligible),
        maximum_artifact_rollout = maximum_artifact_rollout,
        explicit_maximum_artifact_rollout = explicit_maximum_artifact_rollout
    )
}

workflowCapabilityIdentity <- function(
    omic_type,
    workflow_id,
    workflow_type,
    workflow_slug,
    input_format,
    data_level,
    acquisition_mode
) {
    list(
        omic_type = omic_type,
        workflow_id = workflow_id,
        workflow_type = workflow_type,
        workflow_slug = workflow_slug,
        input_format = input_format,
        data_level = data_level,
        acquisition_mode = acquisition_mode
    )
}

.WORKFLOW_CAPABILITY_CATALOGUE <- list(
    newWorkflowCapability(
        "proteomics.diann.peptide.dia.v1",
        workflowCapabilityIdentity(
            "proteomics", "proteomics.gui", "DIA", "prot_dia",
            "diann", "peptide", "dia"
        ),
        artifact_eligible = TRUE,
        auto_eligible = FALSE,
        maximum_artifact_rollout = "dual_write"
    ),
    newWorkflowCapability(
        "proteomics.spectronaut.protein.lfq.v1",
        workflowCapabilityIdentity(
            "proteomics", "proteomics.gui", "LFQ", "prot_spectronaut_lfq",
            "spectronaut", "protein", "dia"
        )
    ),
    newWorkflowCapability(
        "proteomics.spectronaut.peptide.lfq.v1",
        workflowCapabilityIdentity(
            "proteomics", "proteomics.gui", "LFQ", "prot_spectronaut_lfq",
            "spectronaut", "peptide", "dia"
        )
    ),
    newWorkflowCapability(
        "proteomics.fragpipe.protein.lfq.v1",
        workflowCapabilityIdentity(
            "proteomics", "proteomics.gui", "LFQ", "prot_lfq_fragpipe",
            "fragpipe", "protein", "not_recorded"
        )
    ),
    newWorkflowCapability(
        "proteomics.maxquant.protein.lfq.v1",
        workflowCapabilityIdentity(
            "proteomics", "proteomics.gui", "LFQ", "prot_lfq",
            "maxquant", "protein", "not_recorded"
        )
    ),
    newWorkflowCapability(
        "proteomics.pd_tmt.protein.tmt.v1",
        workflowCapabilityIdentity(
            "proteomics", "proteomics.gui", "TMT", "prot_tmt",
            "pd_tmt", "protein", "not_recorded"
        )
    ),
    newWorkflowCapability(
        "metabolomics.msdial.metabolite.standard.v1",
        workflowCapabilityIdentity(
            "metabolomics", "metabolomics.gui", "metabolomics_standard",
            "metab_standard", "msdial", "metabolite", "not_recorded"
        )
    ),
    newWorkflowCapability(
        "metabolomics.custom.metabolite.standard.v1",
        workflowCapabilityIdentity(
            "metabolomics", "metabolomics.gui", "metabolomics_standard",
            "metab_standard", "custom", "metabolite", "not_recorded"
        )
    ),
    newWorkflowCapability(
        "lipidomics.msdial.lipid.standard.v1",
        workflowCapabilityIdentity(
            "lipidomics", "lipidomics.gui", "lipidomics_standard",
            "lipid_standard", "msdial", "lipid", "not_recorded"
        )
    ),
    newWorkflowCapability(
        "lipidomics.lipidsearch.lipid.standard.v1",
        workflowCapabilityIdentity(
            "lipidomics", "lipidomics.gui", "lipidomics_standard",
            "lipid_standard", "lipidsearch", "lipid", "not_recorded"
        )
    ),
    newWorkflowCapability(
        "lipidomics.custom.lipid.standard.v1",
        workflowCapabilityIdentity(
            "lipidomics", "lipidomics.gui", "lipidomics_standard",
            "lipid_standard", "custom", "lipid", "not_recorded"
        )
    )
)
names(.WORKFLOW_CAPABILITY_CATALOGUE) <- vapply(
    .WORKFLOW_CAPABILITY_CATALOGUE,
    `[[`,
    character(1),
    "capability_id"
)
.WORKFLOW_CAPABILITY_INDEX <- stats::setNames(
    seq_along(.WORKFLOW_CAPABILITY_CATALOGUE),
    vapply(
        .WORKFLOW_CAPABILITY_CATALOGUE,
        \(capability) workflowCapabilityKey(capability$identity),
        character(1)
    )
)

workflowCapabilityCatalogue <- function() {
    .WORKFLOW_CAPABILITY_CATALOGUE
}

mergeWorkflowDescriptorCapabilities <- function(
    capabilities = workflowCapabilityCatalogue(),
    descriptor_catalogue = artifactWorkflowDescriptorCatalogue()
) {
    descriptor_capabilities <- artifactDescriptorCapabilities(
        validateArtifactDescriptorCatalogue(descriptor_catalogue)
    )
    if (length(descriptor_capabilities) == 0L) return(capabilities)
    descriptor_keys <- vapply(
        descriptor_capabilities,
        \(capability) workflowCapabilityKey(capability$identity),
        character(1)
    )
    base_keys <- vapply(
        capabilities,
        \(capability) workflowCapabilityKey(capability$identity),
        character(1)
    )
    c(capabilities[!base_keys %in% descriptor_keys], descriptor_capabilities)
}

findWorkflowCapability <- function(identity, capabilities = NULL) {
    key <- workflowCapabilityKey(identity)
    if (is.null(capabilities)) {
        index <- unname(.WORKFLOW_CAPABILITY_INDEX[key])
        if (is.na(index)) return(NULL)
        return(.WORKFLOW_CAPABILITY_CATALOGUE[[index]])
    }
    keys <- vapply(
        capabilities,
        \(capability) workflowCapabilityKey(capability$identity),
        character(1)
    )
    matches <- which(keys == key)
    if (length(matches) > 1L) {
        workflowCapabilityAbort(
            "workflow capability catalogue contains a duplicate exact key",
            "multischolar_duplicate_workflow_capability"
        )
    }
    if (length(matches) == 0L) NULL else capabilities[[matches]]
}

resolveImportedWorkflowIdentity <- function(
    static_identity,
    workflow_type,
    input_format,
    data_level,
    capabilities = NULL
) {
    required_static <- c("project_id", "omic_type", "omic_label", "workflow_id")
    if (!all(required_static %in% names(static_identity)) ||
        !all(vapply(
            static_identity[required_static],
            workflowCapabilityScalarString,
            logical(1)
        )) ||
        !all(vapply(
            list(workflow_type, input_format, data_level),
            workflowCapabilityScalarString,
            logical(1)
        ))) {
        workflowCapabilityAbort(
            "imported workflow identity contains missing or invalid fields",
            "multischolar_invalid_workflow_identity"
        )
    }
    catalogue <- capabilities
    if (is.null(catalogue)) catalogue <- .WORKFLOW_CAPABILITY_CATALOGUE
    matches <- vapply(catalogue, \(capability) {
        candidate <- capability$identity
        identical(candidate$omic_type, static_identity$omic_type) &&
            identical(candidate$workflow_id, static_identity$workflow_id) &&
            identical(candidate$workflow_type, workflow_type) &&
            identical(candidate$input_format, input_format) &&
            identical(candidate$data_level, data_level)
    }, logical(1))
    if (sum(matches) == 1L) {
        return(c(
            static_identity,
            catalogue[[which(matches)]]$identity[
                setdiff(.WORKFLOW_CAPABILITY_KEY_FIELDS, names(static_identity))
            ]
        ))
    }
    if (sum(matches) > 1L) {
        workflowCapabilityAbort(
            "imported workflow identity is ambiguous in the capability catalogue",
            "multischolar_ambiguous_workflow_capability"
        )
    }
    slug_values <- c(static_identity$omic_type, workflow_type, input_format, data_level)
    slug_components <- tolower(gsub("[^A-Za-z0-9]+", "_", slug_values))
    slug <- paste(c("unregistered", slug_components), collapse = "_")
    c(
        static_identity,
        list(
            workflow_type = workflow_type,
            workflow_slug = slug,
            input_format = input_format,
            data_level = data_level,
            acquisition_mode = "not_recorded"
        )
    )
}

normalizeWorkflowProjectState <- function(project_state) {
    if (is.character(project_state)) project_state <- list(status = project_state)
    if (!is.list(project_state) ||
        !workflowCapabilityScalarString(project_state$status) ||
        !project_state$status %in% .WORKFLOW_PROJECT_STATES) {
        workflowCapabilityAbort(
            "workflow project state is invalid",
            "multischolar_invalid_project_state"
        )
    }
    project_state
}

workflowEffectiveRollout <- function(requested_rollout, maximum_rollout, forced) {
    if (is.null(maximum_rollout) ||
        !maximum_rollout %in% .WORKFLOW_ARTIFACT_ROLLOUTS) {
        workflowCapabilityAbort(
            "workflow capability has no certified artifact rollout",
            "multischolar_artifact_not_certified"
        )
    }
    requested_rank <- match(requested_rollout, .WORKFLOW_ARTIFACT_ROLLOUTS)
    maximum_rank <- match(maximum_rollout, .WORKFLOW_ARTIFACT_ROLLOUTS)
    if (isTRUE(forced) && requested_rank > maximum_rank) {
        workflowCapabilityAbort(
            "requested artifact rollout exceeds the certified maximum",
            "multischolar_artifact_rollout_not_certified"
        )
    }
    .WORKFLOW_ARTIFACT_ROLLOUTS[[min(requested_rank, maximum_rank)]]
}

resolveWorkflowBackend <- function(
    identity,
    requested_backend = "auto",
    requested_rollout = "dual_write",
    project_state = "new",
    migration_requested = FALSE,
    capabilities = NULL
) {
    if (!workflowCapabilityScalarString(requested_backend) ||
        !requested_backend %in% .WORKFLOW_BACKENDS) {
        workflowCapabilityAbort(
            "requested workflow backend is invalid",
            "multischolar_invalid_backend_request"
        )
    }
    if (!workflowCapabilityScalarString(requested_rollout) ||
        !requested_rollout %in% .WORKFLOW_ARTIFACT_ROLLOUTS) {
        workflowCapabilityAbort(
            "requested artifact rollout is invalid",
            "multischolar_invalid_rollout_request"
        )
    }
    if (!is.logical(migration_requested) || length(migration_requested) != 1L ||
        is.na(migration_requested)) {
        workflowCapabilityAbort(
            "migration request must be one logical value",
            "multischolar_invalid_migration_request"
        )
    }
    state <- normalizeWorkflowProjectState(project_state)
    capability <- findWorkflowCapability(identity, capabilities)
    capability_id <- if (is.null(capability)) NA_character_ else capability$capability_id
    capability_version <- if (is.null(capability)) {
        NA_character_
    } else {
        capability$capability_version
    }
    result <- \(effective_backend, effective_rollout, reason_code) {
        structure(
            list(
                requested_backend = requested_backend,
                effective_backend = effective_backend,
                effective_rollout = effective_rollout,
                capability_id = capability_id,
                capability_version = capability_version,
                reason_code = reason_code,
                project_state = state$status
            ),
            class = c("MultiScholaRBackendResolution", "list")
        )
    }
    if (identical(requested_backend, "memory")) {
        reason <- if (startsWith(state$status, "artifact_")) {
            "requested_memory_with_artifact_state"
        } else {
            "requested_memory"
        }
        return(result("memory", "none", reason))
    }
    artifact_project <- startsWith(state$status, "artifact_")
    if (identical(requested_backend, "artifact")) {
        if (is.null(capability)) {
            workflowCapabilityAbort(
                "forced artifact mode has no exact workflow capability",
                "multischolar_unknown_workflow_capability"
            )
        }
        if (!isTRUE(capability$artifact_eligible)) {
            workflowCapabilityAbort(
                "forced artifact mode is not certified for this workflow",
                "multischolar_artifact_not_certified"
            )
        }
    }
    if (artifact_project && identical(state$status, "artifact_corrupt")) {
        workflowCapabilityAbort(
            "artifact project state is corrupt and cannot be opened for writing",
            "multischolar_corrupt_artifact_project"
        )
    }
    if (artifact_project && identical(state$status, "artifact_future_schema")) {
        workflowCapabilityAbort(
            "artifact project uses a future unsupported schema",
            "multischolar_future_artifact_schema"
        )
    }
    if (identical(requested_backend, "auto") && artifact_project) {
        if (is.null(capability) || !isTRUE(capability$artifact_eligible)) {
            workflowCapabilityAbort(
                "existing artifact project has no available exact capability",
                "multischolar_artifact_capability_unavailable"
            )
        }
        rollout <- workflowEffectiveRollout(
            requested_rollout,
            capability$maximum_artifact_rollout,
            forced = FALSE
        )
        return(result("artifact", rollout, "auto_existing_artifact"))
    }
    if (identical(state$status, "legacy_memory") && !isTRUE(migration_requested)) {
        if (identical(requested_backend, "artifact")) {
            workflowCapabilityAbort(
                "legacy memory project requires an explicit migration request",
                "multischolar_artifact_migration_required"
            )
        }
        return(result("memory", "none", "auto_preserve_legacy_memory"))
    }
    if (identical(requested_backend, "auto") &&
        (is.null(capability) || !isTRUE(capability$auto_eligible))) {
        reason <- if (is.null(capability)) {
            "auto_unknown_capability_memory"
        } else {
            "auto_capability_not_promoted"
        }
        return(result("memory", "none", reason))
    }
    maximum_rollout <- if (identical(requested_backend, "artifact")) {
        capability$explicit_maximum_artifact_rollout
    } else {
        capability$maximum_artifact_rollout
    }
    rollout <- workflowEffectiveRollout(
        requested_rollout,
        maximum_rollout,
        forced = identical(requested_backend, "artifact")
    )
    reason <- if (identical(requested_backend, "artifact")) {
        "artifact_forced"
    } else if (identical(state$status, "legacy_memory")) {
        "auto_migration_requested"
    } else {
        "auto_promoted_new_project"
    }
    result("artifact", rollout, reason)
}
