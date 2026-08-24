# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

#' Define one exact artifact payload eviction contract
#' @param descriptor Exact immutable workflow descriptor.
#' @param owner_label Human-readable workflow owner.
#' @param owner_condition_name Readiness field proving exact ownership.
#' @param source_fields Workflow fields eligible for release.
#' @param inventory_digest_fn Consumer inventory digest provider.
#' @param compatibility_strategy Rollback reconstruction strategy identifier.
#' @param abort_fn Typed abort function supplied by the workflow adapter.
#' @param invalid_readiness_class Condition class for malformed readiness.
#' @param incomplete_eviction_class Condition class for partial field release.
#' @return A validated payload eviction contract.
#' @noRd
newArtifactPayloadEvictionContract <- function(
    descriptor,
    owner_label,
    owner_condition_name,
    source_fields,
    inventory_digest_fn,
    compatibility_strategy,
    abort_fn,
    invalid_readiness_class,
    incomplete_eviction_class
) {
    descriptor <- validateArtifactWorkflowDescriptor(descriptor)
    valid <- workflowCapabilityScalarString(owner_label) &&
        workflowCapabilityScalarString(owner_condition_name) &&
        is.character(source_fields) && length(source_fields) > 0L &&
        !anyDuplicated(source_fields) && all(nzchar(source_fields)) &&
        is.function(inventory_digest_fn) &&
        workflowCapabilityScalarString(compatibility_strategy) &&
        is.function(abort_fn) &&
        workflowCapabilityScalarString(invalid_readiness_class) &&
        workflowCapabilityScalarString(incomplete_eviction_class)
    if (!isTRUE(valid)) {
        artifactStagePersistenceAbort(
            "artifact payload eviction contract is malformed",
            "multischolar_invalid_artifact_payload_eviction_contract"
        )
    }
    structure(
        list(
            descriptor = descriptor,
            owner_label = owner_label,
            owner_condition_name = owner_condition_name,
            source_fields = source_fields,
            inventory_digest_fn = inventory_digest_fn,
            compatibility_strategy = compatibility_strategy,
            abort_fn = abort_fn,
            invalid_readiness_class = invalid_readiness_class,
            incomplete_eviction_class = incomplete_eviction_class
        ),
        class = "ArtifactPayloadEvictionContract"
    )
}

#' Validate an artifact payload eviction contract
#' @param contract Candidate payload eviction contract.
#' @return The validated contract.
#' @noRd
validateArtifactPayloadEvictionContract <- function(contract) {
    expected <- c(
        "descriptor", "owner_label", "owner_condition_name", "source_fields",
        "inventory_digest_fn", "compatibility_strategy", "abort_fn",
        "invalid_readiness_class", "incomplete_eviction_class"
    )
    valid <- inherits(contract, "ArtifactPayloadEvictionContract") &&
        is.list(contract) && identical(names(contract), expected) &&
        is.function(contract$inventory_digest_fn) && is.function(contract$abort_fn)
    if (!isTRUE(valid)) {
        artifactStagePersistenceAbort(
            "artifact payload eviction contract is invalid",
            "multischolar_invalid_artifact_payload_eviction_contract"
        )
    }
    validateArtifactWorkflowDescriptor(contract$descriptor)
    contract
}

#' Build payload-free identity proofs for immutable artifact refs
#' @param refs Named immutable artifact references.
#' @return Named payload-free reference proofs.
#' @noRd
artifactPayloadReadthroughRefProof <- function(refs) {
    lapply(refs, \(ref) {
        ref <- artifactStoreNormalizeRef(ref)
        list(
            artifact_id = ref$artifact_id,
            generation_id = ref$logical_key$generation_id,
            semantic_digest = ref$hash_policy$semantic$digest,
            byte_digest = ref$hash_policy$byte$digest
        )
    })
}

#' Resolve the effective rollout for an artifact workflow
#' @param workflow_data Mutable workflow state.
#' @param rollout_fn Effective rollout provider.
#' @return The effective rollout, or `NULL` for an unbound workflow.
#' @noRd
artifactPayloadEvictionRollout <- function(workflow_data, rollout_fn) {
    context <- workflow_data$workflow_context
    if (!inherits(context, "WorkflowContext") || !context$isBound()) return(NULL)
    rollout_fn(context)
}

#' Extract required read-through proof flags
#' @param proof Candidate read-through proof.
#' @return A named logical vector.
#' @noRd
artifactPayloadEvictionProofFlags <- function(proof) {
    names <- c(
        "payload_validated", "registry_ready", "current_pointer_verified",
        "readthrough_completed", "annotation_completed"
    )
    stats::setNames(
        vapply(names, \(name) isTRUE(proof[[name]]), logical(1)),
        names
    )
}

#' Validate a rollback checkpoint against read-through proof
#' @param checkpoint Candidate compatibility checkpoint.
#' @param proof Validated read-through proof.
#' @param contract Validated payload eviction contract.
#' @return A scalar logical.
#' @noRd
artifactPayloadEvictionCheckpointReady <- function(
    checkpoint,
    proof,
    contract
) {
    contract <- validateArtifactPayloadEvictionContract(contract)
    descriptor_ready <- is.null(proof$descriptor_id) || (
        identical(checkpoint$descriptor_id, proof$descriptor_id) &&
        identical(proof$descriptor_id, contract$descriptor$descriptor_id)
    )
    isTRUE(descriptor_ready) && is.list(checkpoint) && isTRUE(checkpoint$verified) &&
        identical(checkpoint$strategy, contract$compatibility_strategy) &&
        identical(checkpoint$project_id, proof$project_id) &&
        identical(checkpoint$workflow_id, proof$workflow_id) &&
        identical(checkpoint$design_run_id, proof$design_run_id) &&
        identical(checkpoint$state_generation_id, proof$state_generation_id) &&
        identical(
            checkpoint$consumer_inventory_digest,
            contract$inventory_digest_fn()
        )
}

#' Test whether the state hydration cache is bounded for eviction
#' @param workflow_state Memory or artifact workflow state manager.
#' @return A scalar logical.
#' @noRd
artifactPayloadEvictionCacheReady <- function(workflow_state) {
    if (!inherits(workflow_state, "ArtifactWorkflowState")) return(TRUE)
    info <- workflow_state$getCacheInfo()
    is.list(info) && is.numeric(info$entries) && info$entries <= 1L
}

#' Evaluate descriptor-bound payload eviction readiness
#' @param workflow_data Mutable workflow state.
#' @param contract Validated payload eviction contract.
#' @param rollout_fn Effective rollout provider.
#' @return A named logical readiness vector.
#' @noRd
artifactPayloadEvictionReadiness <- function(
    workflow_data,
    contract,
    rollout_fn = \(context) context$getStorageDecision()$effective_rollout
) {
    contract <- validateArtifactPayloadEvictionContract(contract)
    proof <- workflow_data$artifact_readthrough_proof
    checkpoint <- workflow_data$artifact_compatibility_checkpoint
    source_ready <- vapply(
        contract$source_fields,
        \(name) !is.null(workflow_data[[name]]),
        logical(1)
    )
    owner <- artifactStageCoordinatorOwned(workflow_data, contract$descriptor)
    names(owner) <- contract$owner_condition_name
    c(
        owner,
        evict_rollout = identical(
            artifactPayloadEvictionRollout(workflow_data, rollout_fn),
            "evict"
        ),
        source_payload_available = all(source_ready),
        artifactPayloadEvictionProofFlags(proof),
        compatibility_checkpoint_verified =
            artifactPayloadEvictionCheckpointReady(checkpoint, proof, contract),
        consumer_inventory_verified = identical(
            proof$consumer_inventory_digest,
            contract$inventory_digest_fn()
        ),
        single_cache_entry = artifactPayloadEvictionCacheReady(
            workflow_data$state_manager
        )
    )
}

#' Validate a descriptor-bound readiness vector
#' @param readiness Candidate readiness vector.
#' @param contract Validated payload eviction contract.
#' @return The validated readiness vector.
#' @noRd
validateArtifactPayloadEvictionReadiness <- function(readiness, contract) {
    contract <- validateArtifactPayloadEvictionContract(contract)
    required <- c(
        contract$owner_condition_name, "evict_rollout", "source_payload_available",
        "payload_validated", "registry_ready", "current_pointer_verified",
        "readthrough_completed", "annotation_completed",
        "compatibility_checkpoint_verified", "consumer_inventory_verified",
        "single_cache_entry"
    )
    valid <- is.logical(readiness) && identical(names(readiness), required) &&
        !anyNA(readiness)
    if (!isTRUE(valid)) {
        contract$abort_fn(
            sprintf("%s eviction readiness is malformed", contract$owner_label),
            contract$invalid_readiness_class
        )
    }
    readiness
}

#' Measure source fields eligible for release
#' @param workflow_data Mutable workflow state.
#' @param source_fields Exact source fields.
#' @return Named object sizes in bytes.
#' @noRd
artifactPayloadSourceFieldBytes <- function(workflow_data, source_fields) {
    vapply(
        source_fields,
        \(name) as.numeric(utils::object.size(workflow_data[[name]])),
        numeric(1)
    )
}

#' Assign one workflow field
#' @param workflow_data Mutable workflow state.
#' @param name Field name.
#' @param value Replacement value.
#' @return `TRUE`, invisibly.
#' @noRd
setArtifactWorkflowField <- function(workflow_data, name, value) {
    workflow_data[[name]] <- value
    invisible(TRUE)
}

#' Clear source payload fields as one guarded operation
#' @param workflow_data Mutable workflow state.
#' @param contract Validated payload eviction contract.
#' @param clear_fn Injectable workflow field writer.
#' @param release_cache_fn Injectable hydration-cache release function.
#' @return `TRUE`, invisibly.
#' @noRd
clearArtifactWorkflowPayloads <- function(
    workflow_data,
    contract,
    clear_fn = setArtifactWorkflowField,
    release_cache_fn = workflowStateReleaseHydrationCache
) {
    contract <- validateArtifactPayloadEvictionContract(contract)
    original <- lapply(
        contract$source_fields,
        \(name) workflow_data[[name]]
    )
    names(original) <- contract$source_fields
    committed <- FALSE
    on.exit({
        if (!committed) {
            for (name in contract$source_fields) {
                setArtifactWorkflowField(workflow_data, name, original[[name]])
            }
        }
    }, add = TRUE)
    for (name in contract$source_fields) clear_fn(workflow_data, name, NULL)
    release_cache_fn(workflow_data$state_manager)
    cleared <- vapply(
        contract$source_fields,
        \(name) is.null(workflow_data[[name]]),
        logical(1)
    )
    if (!all(cleared)) {
        contract$abort_fn(
            sprintf(
                "%s payload eviction did not clear every source",
                contract$owner_label
            ),
            contract$incomplete_eviction_class
        )
    }
    committed <- TRUE
    invisible(TRUE)
}

#' Evict descriptor-bound source payload fields
#' @param workflow_data Mutable workflow state.
#' @param contract Validated payload eviction contract.
#' @param rollout_fn Effective rollout provider.
#' @param readiness_fn Optional readiness provider.
#' @param clear_fn Injectable workflow field writer.
#' @param release_cache_fn Injectable hydration-cache release function.
#' @param record_fn Stage result recorder.
#' @return A payload-free eviction result.
#' @noRd
evictArtifactWorkflowPayloads <- function(
    workflow_data,
    contract,
    rollout_fn = \(context) context$getStorageDecision()$effective_rollout,
    readiness_fn = NULL,
    clear_fn = setArtifactWorkflowField,
    release_cache_fn = workflowStateReleaseHydrationCache,
    record_fn = recordArtifactStageResult
) {
    contract <- validateArtifactPayloadEvictionContract(contract)
    previous <- workflow_data$artifact_stage_results$eviction
    if (is.list(previous) && isTRUE(previous$evicted)) return(previous)
    readiness <- if (is.null(readiness_fn)) {
        artifactPayloadEvictionReadiness(workflow_data, contract, rollout_fn)
    } else {
        readiness_fn(workflow_data, rollout_fn)
    }
    readiness <- validateArtifactPayloadEvictionReadiness(readiness, contract)
    failed <- names(readiness)[!readiness]
    if (length(failed) > 0L) {
        return(list(
            enabled = isTRUE(readiness[[contract$owner_condition_name]]) &&
                isTRUE(readiness[["evict_rollout"]]),
            ok = TRUE,
            evicted = FALSE,
            reason = "eviction_prerequisites_incomplete",
            failed_prerequisites = failed
        ))
    }
    source_field_bytes <- artifactPayloadSourceFieldBytes(
        workflow_data,
        contract$source_fields
    )
    clearArtifactWorkflowPayloads(
        workflow_data,
        contract,
        clear_fn,
        release_cache_fn
    )
    result <- list(
        enabled = TRUE,
        ok = TRUE,
        evicted = TRUE,
        reason = "artifact_payloads_evicted",
        descriptor_id = contract$descriptor$descriptor_id,
        evicted_fields = contract$source_fields,
        released_source_field_bytes = source_field_bytes,
        released_source_bytes_upper_bound = sum(source_field_bytes),
        state_generation_id =
            workflow_data$artifact_readthrough_proof$state_generation_id,
        compatibility_strategy = contract$compatibility_strategy
    )
    record_fn(workflow_data, "eviction", result)
}

#' Resolve one explicitly requested evicted source table
#' @param workflow_data Mutable workflow state.
#' @param name Exact source field name.
#' @param contract Validated payload eviction contract.
#' @param ref_fn Function resolving the source field artifact reference.
#' @param read_fn Function independently hydrating one artifact reference.
#' @return The retained table, reconstructed table, or `NULL`.
#' @noRd
resolveArtifactWorkflowTable <- function(
    workflow_data,
    name,
    contract,
    ref_fn,
    read_fn
) {
    contract <- validateArtifactPayloadEvictionContract(contract)
    if (!name %in% contract$source_fields) {
        contract$abort_fn(
            sprintf(
                "%s payload resolver received an unsupported field",
                contract$owner_label
            ),
            contract$invalid_readiness_class
        )
    }
    if (!is.null(workflow_data[[name]])) return(workflow_data[[name]])
    if (!artifactStageCoordinatorOwned(workflow_data, contract$descriptor)) {
        return(NULL)
    }
    ref <- ref_fn(workflow_data, name)
    if (!is.list(ref)) return(NULL)
    identity <- workflow_data$workflow_context$getIdentity()
    store <- newArtifactStore(
        workflow_data$workflow_context$getPaths(),
        identity$project_id
    )
    read_fn(store, ref)
}
