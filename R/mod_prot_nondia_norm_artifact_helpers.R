# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

.PROT_NONDIA_NORM_AUDIT_SCHEMA <-
    "multischolar.prot_nondia.normalization_audit"
.PROT_NONDIA_NORM_AUDIT_VERSION <- 1L
.PROT_NONDIA_NORM_REF_SCHEMA <-
    "multischolar.prot_nondia.normalization_state_ref"
.PROT_NONDIA_NORM_REF_VERSION <- 1L

#' Test whether coordinator-owned proteomics workflow data is available
#' @param workflow_data Candidate mutable workflow state.
#' @return A scalar logical.
#' @noRd
protNormWorkflowDataAvailable <- function(workflow_data) {
    if (is.null(workflow_data)) return(FALSE)
    manager <- tryCatch(workflow_data$state_manager, error = \(...) NULL)
    !is.null(manager)
}

#' Resolve the independent artifact switch for one normalization stage
#' @param workflow_data Candidate mutable workflow state.
#' @param stage_id Optional normalization stage identifier.
#' @return One of `"enabled"` or `"disabled"`.
#' @noRd
protNormArtifactMode <- function(workflow_data, stage_id = NULL) {
    if (!is.null(protNonDiaNormDescriptor(workflow_data))) {
        return(protNonDiaNormArtifactMode(stage_id))
    }
    protDiaNormArtifactMode(stage_id)
}

#' Resolve the independent non-DIA normalization artifact switch
#' @param stage_id Optional normalization stage identifier.
#' @return One of `"enabled"` or `"disabled"`.
#' @noRd
protNonDiaNormArtifactMode <- function(stage_id = NULL) {
    mode <- getOption(
        "multischolar.prot_nondia.normalization_artifacts",
        "enabled"
    )
    if (!is.null(stage_id)) {
        stage_id <- match.arg(stage_id, .PROT_DIA_NORM_STAGES)
        mode <- getOption(
            paste0(
                "multischolar.prot_nondia.normalization_artifacts.",
                stage_id
            ),
            mode
        )
    }
    match.arg(mode, c("enabled", "disabled"))
}

#' Resolve an exact non-DIA descriptor from workflow state
#' @param workflow_data Mutable proteomics workflow state.
#' @return The exact supported descriptor, or `NULL` when ineligible.
#' @noRd
protNonDiaNormDescriptor <- function(workflow_data) {
    descriptor <- tryCatch(
        protArtifactBoundDescriptor(workflow_data),
        error = \(...) NULL
    )
    if (is.null(descriptor) ||
        !descriptor$descriptor_id %in% names(protNonDiaReadthroughDescriptors())) {
        return(NULL)
    }
    protNonDiaReadthroughDescriptor(descriptor$descriptor_id)
}

#' Test exact non-DIA normalization artifact eligibility
#' @param workflow_data Mutable proteomics workflow state.
#' @param state_manager Current workflow state manager.
#' @param before Parent scientific S4 state.
#' @param after Candidate scientific S4 state.
#' @param stage_id Normalization stage identifier.
#' @return A scalar logical.
#' @noRd
protNonDiaNormArtifactEligible <- function(
    workflow_data,
    state_manager,
    before,
    after,
    stage_id
) {
    if (!protDiaNormTransition(before, after, stage_id) ||
        identical(protNonDiaNormArtifactMode(stage_id), "disabled")) {
        return(FALSE)
    }
    descriptor <- protNonDiaNormDescriptor(workflow_data)
    if (is.null(descriptor)) return(FALSE)
    workflow_type <- tryCatch(
        workflowStateType(state_manager),
        error = \(...) NULL
    )
    identical(workflow_type, descriptor$identity$workflow_type) &&
        artifactStageCoordinatorOwned(workflow_data, descriptor)
}

#' Construct or reuse the exact non-DIA normalization state manager
#' @param workflow_data Mutable proteomics workflow state.
#' @param state_manager Current memory or artifact state manager.
#' @param descriptor Exact supported workflow descriptor.
#' @return A manager resource record with an ownership flag.
#' @noRd
protNonDiaNormArtifactManager <- function(
    workflow_data,
    state_manager,
    descriptor
) {
    if (inherits(state_manager, "ArtifactWorkflowState")) {
        return(list(manager = state_manager, owned = FALSE))
    }
    list(
        manager = newWorkflowState(
            workflow_context = workflow_data$workflow_context,
            workflow_descriptor = descriptor,
            descriptor_catalogue = artifactWorkflowDescriptorCatalogue(),
            codec_catalogue = artifactS4CodecCatalogue()
        ),
        owned = TRUE
    )
}

#' Validate the active non-DIA normalization parent state
#' @param manager Descriptor-bound artifact state manager.
#' @param before Expected parent S4 state.
#' @return The validated parent, invisibly.
#' @noRd
protNonDiaNormValidateParent <- function(manager, before) {
    parent <- manager$getState()
    if (!identical(parent, before)) {
        protNonDiaArtifactAbort(
            "non-DIA normalization parent differs from the active scientific state",
            "multischolar_prot_nondia_normalization_parent_mismatch"
        )
    }
    invisible(parent)
}

#' Resolve the current non-DIA normalization audit record
#' @param object Protein S4 state.
#' @return The current audit record, or `NULL` when no audit exists.
#' @noRd
protNonDiaNormCurrentRecord <- function(object) {
    audit <- protDiaNormArgs(object)$normalization_artifact_audit
    if (!is.list(audit) || !is.list(audit$records) ||
        length(audit$records) == 0L) {
        return(NULL)
    }
    record <- tail(audit$records, 1L)[[1L]]
    if (!identical(record$record_id, audit$current_record_id)) {
        protNonDiaArtifactAbort(
            "non-DIA normalization audit head differs from its latest record",
            "multischolar_invalid_prot_nondia_normalization_audit"
        )
    }
    record
}

#' Emit one tuple-bound non-DIA normalization audit record
#' @param before Parent protein S4 state.
#' @param after Candidate protein S4 state.
#' @param descriptor Exact supported workflow descriptor.
#' @param stage_id Normalization stage identifier.
#' @param parent_state Parent state name.
#' @param current_state Candidate state name.
#' @param parent_generation_id Parent generation identifier.
#' @param parameters Exact scientific parameters.
#' @param status Applied or skipped status.
#' @param transformation_type Scientific transformation type.
#' @param software Software provenance.
#' @param now Audit timestamp.
#' @return The candidate S4 state with appended audit metadata.
#' @noRd
protNonDiaNormEmitAudit <- function(
    before,
    after,
    descriptor,
    stage_id,
    parent_state,
    current_state,
    parent_generation_id,
    parameters,
    status,
    transformation_type,
    software,
    now
) {
    if (!methods::is(after, "ProteinQuantitativeData")) return(after)
    prior <- protNonDiaNormCurrentRecord(before)
    parent_record_id <- if (is.null(prior)) {
        protDiaNormInheritedRecordId(
            before,
            "protein_qc_audit",
            "legacy_untracked_protein_qc"
        )
    } else {
        prior$record_id
    }
    substantive <- list(
        schema = .PROT_NONDIA_NORM_AUDIT_SCHEMA,
        schema_version = .PROT_NONDIA_NORM_AUDIT_VERSION,
        capability_id = descriptor$descriptor_id,
        capability_version = descriptor$descriptor_version,
        stage_id = stage_id,
        status = status,
        transformation_type = transformation_type,
        parent_state = parent_state,
        current_state = current_state,
        parent_generation_id = parent_generation_id,
        parent_record_id = parent_record_id,
        peptide_qc_record_id = protDiaNormInheritedRecordId(
            after,
            "peptide_qc_audit",
            "not_reachable_non_dia"
        ),
        protein_qc_record_id = protDiaNormInheritedRecordId(
            after,
            "protein_qc_audit",
            "not_reachable_non_dia"
        ),
        resolved_parameters = protDiaNormAuditParameters(parameters),
        software = software,
        before_summary = protDiaNormSummary(before),
        after_summary = protDiaNormSummary(after)
    )
    record_id <- paste0(
        "prot-nondia-norm:",
        substr(.peptideQcDigest(substantive), 1L, 24L)
    )
    record <- c(
        substantive,
        list(
            record_id = record_id,
            timestamp_utc = format(
                now,
                "%Y-%m-%dT%H:%M:%OS6Z",
                tz = "UTC"
            ),
            canonical_digest = .peptideQcDigest(substantive)
        )
    )
    previous <- protDiaNormArgs(before)$normalization_artifact_audit$records
    if (!is.list(previous)) previous <- list()
    previous[[length(previous) + 1L]] <- record
    after@args$normalization_artifact_audit <- list(
        schema = .PROT_NONDIA_NORM_AUDIT_SCHEMA,
        schema_version = .PROT_NONDIA_NORM_AUDIT_VERSION,
        records = previous,
        current_record_id = record_id
    )
    after
}

#' Save one exact non-DIA normalization artifact state
#' @param workflow_data Mutable proteomics workflow state.
#' @param state_manager Current memory or artifact state manager.
#' @param before Parent protein S4 state.
#' @param after Candidate protein S4 state.
#' @param stage_id Normalization stage identifier.
#' @param state_name Candidate state name.
#' @param config_object State configuration.
#' @param description State description.
#' @param parameters Exact scientific parameters.
#' @param status Applied or skipped status.
#' @param transformation_type Scientific transformation type.
#' @param now Audit timestamp.
#' @param failure_injector Optional artifact failure injector used by tests.
#' @return A structured artifact persistence result.
#' @noRd
saveProtNonDiaNormArtifactState <- function(
    workflow_data,
    state_manager,
    before,
    after,
    stage_id,
    state_name,
    config_object,
    description,
    parameters,
    status,
    transformation_type,
    now,
    failure_injector = NULL
) {
    if (!protNonDiaNormArtifactEligible(
        workflow_data,
        state_manager,
        before,
        after,
        stage_id
    )) {
        return(list(handled = FALSE, enabled = FALSE, object = after))
    }
    descriptor <- protNonDiaNormDescriptor(workflow_data)
    resources <- protNonDiaNormArtifactManager(
        workflow_data,
        state_manager,
        descriptor
    )
    manager <- resources$manager
    if (isTRUE(resources$owned)) on.exit(manager$close(), add = TRUE)
    protNonDiaNormValidateParent(manager, before)
    parent_state <- manager$getCurrentStateName()
    memory_parent_state <- workflowStateCurrentName(state_manager)
    if (isTRUE(resources$owned) && !identical(parent_state, memory_parent_state)) {
        protNonDiaArtifactAbort(
            "non-DIA normalization managers disagree on the active parent state",
            "multischolar_prot_nondia_normalization_parent_mismatch"
        )
    }
    parent_generation_id <- manager$getCurrentGenerationId()
    after <- protDiaNormAttachParameters(after, parameters, stage_id)
    software <- protDiaNormSoftware(stage_id, parameters)
    after <- protNonDiaNormEmitAudit(
        before,
        after,
        descriptor,
        stage_id,
        parent_state,
        state_name,
        parent_generation_id,
        parameters,
        status,
        transformation_type,
        software,
        now
    )
    record <- protNonDiaNormCurrentRecord(after)
    config_object$audit_record_id <- record$record_id
    artifactStoreInvokeFailure(
        failure_injector,
        "after_audit_creation",
        list(
            stage_id = stage_id,
            state_name = state_name,
            audit_record_id = record$record_id,
            parent_generation_id = parent_generation_id
        )
    )
    result <- manager$commitState(
        state_name = state_name,
        s4_data_object = after,
        config_object = protDiaNormSafeConfig(config_object, software),
        description = description,
        audit_metadata = record,
        failure_injector = failure_injector,
        action_id = artifactOpaqueId("action"),
        expected_parent_generation_id = parent_generation_id
    )
    if (!identical(result$status, "accepted")) {
        protNonDiaArtifactAbort(
            "non-DIA normalization state did not advance its exact parent",
            "multischolar_prot_nondia_normalization_commit_rejected",
            result = result
        )
    }
    if (isTRUE(resources$owned)) {
        memory_error <- tryCatch(
            {
                protDiaNormSaveMemory(
                    state_manager,
                    state_name,
                    after,
                    config_object,
                    description,
                    record
                )
                NULL
            },
            error = identity
        )
        if (inherits(memory_error, "error")) {
            restore_error <- tryCatch(
                {
                    manager$revertToState(parent_state)
                    NULL
                },
                error = identity
            )
            if (inherits(restore_error, "error")) {
                protNonDiaArtifactAbort(
                    paste(
                        "non-DIA normalization compatibility save and",
                        "artifact rollback both failed"
                    ),
                    "multischolar_prot_nondia_normalization_divergent_state",
                    memory_error = memory_error,
                    restore_error = restore_error
                )
            }
            stop(memory_error)
        }
    }
    hydrated <- manager$getState(state_name)
    role <- artifactProteomicsNonDiaCodecRole(
        descriptor$descriptor_id,
        state_name,
        descriptor$descriptor_version
    )
    artifactValidateProteomicsNonDiaProteinState(hydrated, role, state_name)
    if (!identical(hydrated, after) ||
        !isTRUE(methods::validObject(hydrated, test = TRUE))) {
        protNonDiaArtifactAbort(
            "non-DIA normalization artifact hydration changed scientific state",
            "multischolar_inexact_prot_nondia_normalization_hydration"
        )
    }
    output <- list(
        handled = TRUE,
        enabled = TRUE,
        ok = TRUE,
        committed = TRUE,
        object = after,
        capability_id = descriptor$descriptor_id,
        stage_id = stage_id,
        state_name = state_name,
        generation_id = result$generation_id,
        parent_generation_id = parent_generation_id,
        representation = result$representation,
        metrics = result$metrics,
        audit_record_id = record$record_id
    )
    recordArtifactStageResult(workflow_data, stage_id, output)
    output
}

#' Build one payload-free non-DIA normalization state reference
#' @param workflow_data Mutable proteomics workflow state.
#' @param stage_id Normalization stage identifier.
#' @param state_name State name.
#' @return A generation reference, or `NULL` when no artifact was committed.
#' @noRd
protNonDiaNormStateRef <- function(workflow_data, stage_id, state_name) {
    result <- workflow_data$artifact_stage_results[[stage_id]]
    if (!is.list(result) || !isTRUE(result$enabled) ||
        !workflowStateScalarString(result$generation_id) ||
        !workflowStateScalarString(result$capability_id)) {
        return(NULL)
    }
    list(
        schema = .PROT_NONDIA_NORM_REF_SCHEMA,
        schema_version = .PROT_NONDIA_NORM_REF_VERSION,
        capability_id = result$capability_id,
        stage_id = stage_id,
        state_name = state_name,
        generation_id = result$generation_id
    )
}

#' Revert an exact non-DIA normalization state across dual managers
#' @param workflow_data Mutable proteomics workflow state.
#' @param state_name Target state name.
#' @return The exact reverted memory state.
#' @noRd
revertProtNonDiaNormState <- function(workflow_data, state_name) {
    state_manager <- workflow_data$state_manager
    if (inherits(state_manager, "ArtifactWorkflowState")) {
        return(state_manager$revertToState(state_name))
    }
    descriptor <- protNonDiaNormDescriptor(workflow_data)
    if (is.null(descriptor) ||
        identical(protNonDiaNormArtifactMode(), "disabled")) {
        return(state_manager$revertToState(state_name))
    }
    artifact_manager <- protNonDiaNormArtifactManager(
        workflow_data,
        state_manager,
        descriptor
    )$manager
    on.exit(artifact_manager$close(), add = TRUE)
    memory_current <- state_manager$getCurrentStateName()
    artifact_current <- artifact_manager$getCurrentStateName()
    targets_match <- isTRUE(state_manager$hasState(state_name)) &&
        isTRUE(artifact_manager$hasState(state_name)) &&
        identical(memory_current, artifact_current) &&
        identical(state_manager$getState(), artifact_manager$getState()) &&
        identical(
            state_manager$getState(state_name),
            artifact_manager$getState(state_name)
        )
    if (!targets_match) {
        protNonDiaArtifactAbort(
            "non-DIA normalization managers differ before revert",
            "multischolar_prot_nondia_normalization_parent_mismatch"
        )
    }
    reverted_artifact <- artifact_manager$revertToState(state_name)
    memory_error <- tryCatch(
        {
            state_manager$revertToState(state_name)
            NULL
        },
        error = identity
    )
    if (inherits(memory_error, "error")) {
        protDiaPeptideQcRestoreRevert(
            artifact_manager,
            state_manager,
            artifact_current,
            memory_current
        )
        stop(memory_error)
    }
    reverted_memory <- state_manager$getState()
    if (!identical(reverted_artifact, reverted_memory)) {
        protDiaPeptideQcRestoreRevert(
            artifact_manager,
            state_manager,
            artifact_current,
            memory_current
        )
        protNonDiaArtifactAbort(
            "non-DIA normalization revert hydrated different active states",
            "multischolar_inexact_prot_nondia_normalization_hydration"
        )
    }
    reverted_memory
}

#' Dispatch normalization revert by exact proteomics descriptor
#' @param workflow_data Mutable proteomics workflow state.
#' @param state_name Target state name.
#' @return The exact reverted state.
#' @noRd
revertProtNormState <- function(workflow_data, state_name) {
    if (!is.null(protNonDiaNormDescriptor(workflow_data))) {
        return(revertProtNonDiaNormState(workflow_data, state_name))
    }
    revertProtDiaNormState(workflow_data, state_name)
}
