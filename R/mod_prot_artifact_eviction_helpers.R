# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

PROT_DIA_EVICT_FIELDS <- c("data_tbl", "data_cln")
PROT_DIA_EVICTION_PROOF_SCHEMA <- "multischolar.prot_dia_eviction_proof"
PROT_DIA_EVICTION_PROOF_VERSION <- 1L

#' Define the exact DIA-NN payload eviction compatibility contract
#' @return A descriptor-bound artifact payload eviction contract.
#' @noRd
protDiaPayloadEvictionContract <- function() {
    newArtifactPayloadEvictionContract(
        descriptor = artifactDiaWorkflowDescriptor(),
        owner_label = "DIA-NN",
        owner_condition_name = "exact_dia_canary",
        source_fields = PROT_DIA_EVICT_FIELDS,
        inventory_digest_fn = protDiaEvictionInventoryDigest,
        compatibility_strategy = "read_only_artifact_reconstruction",
        abort_fn = protDiaArtifactAbort,
        invalid_readiness_class =
            "multischolar_invalid_prot_dia_eviction_readiness",
        incomplete_eviction_class =
            "multischolar_incomplete_prot_dia_eviction"
    )
}

#' Inventory DIA-NN consumers at the import/design eviction boundary
#'
#' @return A data frame declaring each consumer and its post-eviction source.
#' @noRd
protDiaEvictionConsumerInventory <- function() {
    readers <- c(
        "registerProtImportOutputs",
        "registerProtDesignPreviewOutputs",
        "registerProtDesignBuilderModule",
        "resolveProtEnrichUniprotAnnotations",
        "hydrateArtifactStageInput",
        "collectProtNormExportSessionData",
        "prepareProtSummaryCopyInputs",
        "resumeProtDiaArtifactWorkflow",
        "resolveProtDiaWorkflowTable"
    )
    data.frame(
        category = c(
            "preview", "preview", "design", "annotation", "qc", "session",
            "report", "compatibility", "compatibility"
        ),
        reader = readers,
        source_after_eviction = c(
            "bounded import ref and setup log",
            "retained import metadata and immutable import ref",
            "canonical import ref for explicit design rebuild",
            "retained annotations/sequences and current S4 state",
            "current workflow state generation",
            "current workflow state plus retained metadata",
            "current workflow state, retained metadata, compatibility files",
            "read-only immutable import/design refs and state manifest",
            "one explicitly requested immutable table ref"
        ),
        verified = vapply(readers, \(reader) {
            exists(
                reader,
                envir = environment(protDiaEvictionConsumerInventory),
                mode = "function",
                inherits = TRUE
            )
        }, logical(1)),
        stringsAsFactors = FALSE
    )
}

protDiaEvictionInventoryDigest <- function() {
    artifactSemanticDigest(protDiaEvictionConsumerInventory())
}

protDiaReadthroughRefProof <- function(refs) {
    artifactPayloadReadthroughRefProof(refs)
}

#' Build a small proof after exact DIA-NN artifact hydration
#'
#' @param bundle A validated DIA-NN resume bundle.
#' @return A payload-free read-through proof.
#' @noRd
newProtDiaReadthroughProof <- function(bundle) {
    evidence <- bundle$evidence
    list(
        schema = PROT_DIA_EVICTION_PROOF_SCHEMA,
        schema_version = PROT_DIA_EVICTION_PROOF_VERSION,
        project_id = evidence$identity$project_id,
        workflow_id = evidence$identity$workflow_id,
        import_run_id = evidence$import$run_id,
        design_run_id = evidence$design$run_id,
        state_generation_id = evidence$current_state$generation_id,
        payload_validated = TRUE,
        registry_ready = TRUE,
        current_pointer_verified = TRUE,
        readthrough_completed = TRUE,
        annotation_completed = isTRUE(bundle$annotation_completed),
        readthrough_mode = bundle$readthrough_mode,
        import_refs = protDiaReadthroughRefProof(evidence$import$refs),
        design_refs = protDiaReadthroughRefProof(evidence$design$refs),
        consumer_inventory_digest = protDiaEvictionInventoryDigest(),
        verified_at = artifactRefUtcNow()
    )
}

newProtDiaCompatibilityCheckpoint <- function(proof) {
    list(
        strategy = "read_only_artifact_reconstruction",
        verified = TRUE,
        project_id = proof$project_id,
        workflow_id = proof$workflow_id,
        design_run_id = proof$design_run_id,
        state_generation_id = proof$state_generation_id,
        consumer_inventory_digest = proof$consumer_inventory_digest,
        verified_at = proof$verified_at
    )
}

recordProtDiaReadthroughProof <- function(workflow_data, bundle) {
    proof <- newProtDiaReadthroughProof(bundle)
    checkpoint <- newProtDiaCompatibilityCheckpoint(proof)
    workflow_data$artifact_readthrough_proof <- proof
    workflow_data$artifact_compatibility_checkpoint <- checkpoint
    list(proof = proof, compatibility_checkpoint = checkpoint)
}

protDiaEvictionRollout <- function(workflow_data, rollout_fn) {
    artifactPayloadEvictionRollout(workflow_data, rollout_fn)
}

protDiaEvictionProofFlags <- function(proof) {
    artifactPayloadEvictionProofFlags(proof)
}

protDiaEvictionCheckpointReady <- function(checkpoint, proof) {
    artifactPayloadEvictionCheckpointReady(
        checkpoint,
        proof,
        protDiaPayloadEvictionContract()
    )
}

protDiaEvictionCacheReady <- function(workflow_state) {
    artifactPayloadEvictionCacheReady(workflow_state)
}

#' Evaluate all DIA-NN eviction prerequisites without mutation
#'
#' @param workflow_data The proteomics workflow state bus.
#' @param rollout_fn Function resolving the effective artifact rollout.
#' @return A named logical vector.
#' @noRd
protDiaEvictionReadiness <- function(
    workflow_data,
    rollout_fn = \(context) context$getStorageDecision()$effective_rollout
) {
    artifactPayloadEvictionReadiness(
        workflow_data,
        protDiaPayloadEvictionContract(),
        rollout_fn
    )
}

validateProtDiaEvictionReadiness <- function(readiness) {
    validateArtifactPayloadEvictionReadiness(
        readiness,
        protDiaPayloadEvictionContract()
    )
}

protDiaEvictionSourceFieldBytes <- function(workflow_data) {
    artifactPayloadSourceFieldBytes(
        workflow_data,
        PROT_DIA_EVICT_FIELDS
    )
}

setProtDiaWorkflowField <- function(workflow_data, name, value) {
    setArtifactWorkflowField(workflow_data, name, value)
}

clearProtDiaWorkflowPayloads <- function(
    workflow_data,
    clear_fn,
    release_cache_fn
) {
    clearArtifactWorkflowPayloads(
        workflow_data,
        protDiaPayloadEvictionContract(),
        clear_fn,
        release_cache_fn
    )
}

#' Evict settled DIA-NN import/design reactive payloads
#'
#' @param workflow_data The proteomics workflow state bus.
#' @param rollout_fn Function resolving the effective artifact rollout.
#' @param readiness_fn Function evaluating all release prerequisites.
#' @param clear_fn Function assigning a workflow field, injectable for fault tests.
#' @param release_cache_fn Function releasing the state hydration cache.
#' @param gc_fn Function collecting unreachable R objects after commit.
#' @return A payload-free eviction result.
#' @noRd
evictProtDiaWorkflowPayloads <- function(
    workflow_data,
    rollout_fn = \(context) context$getStorageDecision()$effective_rollout,
    readiness_fn = protDiaEvictionReadiness,
    clear_fn = setProtDiaWorkflowField,
    release_cache_fn = workflowStateReleaseHydrationCache,
    gc_fn = NULL
) {
    result <- evictArtifactWorkflowPayloads(
        workflow_data,
        protDiaPayloadEvictionContract(),
        rollout_fn,
        readiness_fn,
        clear_fn,
        release_cache_fn
    )
    result
}

evictProtDiaWorkflowPayloadsSafely <- function(
    workflow_data,
    rollout_fn = \(context) context$getStorageDecision()$effective_rollout,
    log_warn = logger::log_warn,
    ...
) {
    result <- tryCatch(
        evictProtDiaWorkflowPayloads(
            workflow_data,
            rollout_fn = rollout_fn,
            ...
        ),
        error = \(error) {
            log_warn(paste(
                "DIA-NN payload eviction failed; source values were restored:",
                conditionMessage(error)
            ))
            list(
                enabled = TRUE,
                ok = FALSE,
                evicted = FALSE,
                reason = "artifact_eviction_failed",
                error_class = class(error),
                error_message = conditionMessage(error)
            )
        }
    )
    recordProtDiaArtifactResult(workflow_data, "eviction", result)
}

verifyProtDiaReadthroughForEviction <- function(
    workflow_data,
    experiment_paths,
    experiment_label,
    storage_policy = NULL,
    resource_policy = NULL,
    parity_process_fn = runProtDiaReadthroughParityProcess,
    parity_failure_stage = NULL
) {
    result <- if (protDiaStageParityAvailable(workflow_data) &&
        is.null(parity_failure_stage)) {
        verifyProtDiaBoundedStageParity(workflow_data, resource_policy)
    } else {
        parity_process_fn(protDiaParityWorkerSpec(
            experiment_paths,
            experiment_label,
            storage_policy,
            resource_policy,
            failure_stage = parity_failure_stage
        ))
    }
    workflow_data$artifact_readthrough_refs <- list(
        import = result$import_refs,
        design = result$design_refs
    )
    workflow_data$artifact_readthrough_proof <- result$proof
    workflow_data$artifact_compatibility_checkpoint <-
        result$compatibility_checkpoint
    list(
        enabled = TRUE,
        ok = TRUE,
        resumed = FALSE,
        verified = TRUE,
        reason = "artifact_readthrough_verified_for_eviction",
        project_id = result$project_id,
        import_run_id = result$import_run_id,
        design_run_id = result$design_run_id,
        state_generation_id = result$state_generation_id,
        compatibility_checkpoint = result$compatibility_checkpoint,
        parity_worker_pid = result$worker_pid,
        complete_payload_returned = FALSE,
        parity_reused = isTRUE(result$parity_reused)
    )
}

newProtDiaSettledStateManagerFromContext <- function(
    context,
    resource_policy = NULL,
    state_export = NULL
) {
    if (!is.list(state_export)) {
        protDiaArtifactAbort(
            "DIA-NN settlement lacks a committed state export",
            "multischolar_missing_prot_dia_settled_state"
        )
    }
    identity <- context$getIdentity()
    store <- newArtifactStore(context$getPaths(), identity$project_id)
    bootstrap <- artifactWorkflowStateBootstrapFromExport(store, state_export)
    decision <- context$getStorageDecision()
    descriptor <- if (identical(decision$capability_version, "1.0.0")) {
        artifactDiaWorkflowDescriptorV1()
    } else {
        artifactDiaWorkflowDescriptor()
    }
    newWorkflowState(
        workflow_context = context,
        resource_policy = resource_policy,
        workflow_descriptor = descriptor,
        descriptor_catalogue = if (identical(
            descriptor$descriptor_version,
            artifactDiaWorkflowDescriptor()$descriptor_version
        )) {
            artifactWorkflowDescriptorCatalogue()
        } else {
            NULL
        },
        codec_catalogue = artifactS4CodecCatalogue(),
        settled_bootstrap = bootstrap
    )
}

newProtDiaSettledStateManager <- function(
    workflow_data,
    resource_policy = NULL,
    state_export = NULL
) {
    if (is.null(state_export)) {
        state_export <- workflow_data$artifact_stage_results$design$state_manifest
    }
    newProtDiaSettledStateManagerFromContext(
        workflow_data$workflow_context,
        resource_policy,
        state_export
    )
}

settleProtDiaArtifactWorkflowSafely <- function(
    workflow_data,
    experiment_paths,
    experiment_label,
    storage_policy = NULL,
    resource_policy = NULL,
    rollout_fn = \(context) context$getStorageDecision()$effective_rollout,
    log_warn = logger::log_warn,
    parity_process_fn = runProtDiaReadthroughParityProcess,
    parity_failure_stage = NULL
) {
    rollout <- protDiaEvictionRollout(workflow_data, rollout_fn)
    if (!identical(rollout, "evict")) {
        return(list(
            enabled = FALSE,
            ok = TRUE,
            evicted = FALSE,
            reason = "rollout_below_evict"
        ))
    }
    result <- tryCatch(
        {
            if (is.null(workflow_data$artifact_readthrough_proof)) {
                verification <- verifyProtDiaReadthroughForEviction(
                    workflow_data,
                    experiment_paths,
                    experiment_label,
                    storage_policy,
                    resource_policy,
                    parity_process_fn,
                    parity_failure_stage
                )
            } else {
                verification <- list(verified = TRUE)
            }
            if (!isTRUE(verification$verified)) {
                verification
            } else {
                settled_manager <- newProtDiaSettledStateManager(
                    workflow_data,
                    resource_policy
                )
                manager_installed <- FALSE
                on.exit({
                    if (!manager_installed) {
                        try(settled_manager$close(), silent = TRUE)
                    }
                }, add = TRUE)
                previous_manager <- workflow_data$state_manager
                eviction <- evictProtDiaWorkflowPayloads(
                    workflow_data,
                    rollout_fn = rollout_fn
                )
                workflow_data$state_manager <- settled_manager
                manager_installed <- TRUE
                if (inherits(previous_manager, "ArtifactWorkflowState")) {
                    try(previous_manager$close(), silent = TRUE)
                }
                eviction$parity_worker_pid <- verification$parity_worker_pid
                eviction$complete_payload_returned <-
                    verification$complete_payload_returned
                eviction$parity_reused <- verification$parity_reused
                eviction$state_manager_replaced <- TRUE
                eviction
            }
        },
        error = \(error) {
            log_warn(paste(
                "DIA-NN artifact settlement failed without evicting sources:",
                conditionMessage(error)
            ))
            list(
                enabled = TRUE,
                ok = FALSE,
                evicted = FALSE,
                reason = "artifact_settlement_failed",
                error_class = class(error),
                error_message = conditionMessage(error)
            )
        }
    )
    recordProtDiaArtifactResult(workflow_data, "eviction", result)
}

protDiaWorkflowPayloadRef <- function(workflow_data, name) {
    refs <- workflow_data$artifact_readthrough_refs
    if (identical(name, "data_tbl") && is.list(refs$import$canonical_data)) {
        return(refs$import$canonical_data)
    }
    if (identical(name, "data_cln") && is.list(refs$design$cleaned_data)) {
        return(refs$design$cleaned_data)
    }
    stages <- workflow_data$artifact_stage_results
    if (identical(name, "data_cln") && is.list(stages$design$refs$cleaned_data)) {
        return(stages$design$refs$cleaned_data)
    }
    stages$import$refs$canonical_data
}

protDiaWorkflowPayloadAvailable <- function(workflow_data, name) {
    if (!name %in% PROT_DIA_EVICT_FIELDS) return(!is.null(workflow_data[[name]]))
    if (!is.null(workflow_data[[name]])) return(TRUE)
    is.list(protDiaWorkflowPayloadRef(workflow_data, name))
}

protDiaWorkflowPayloadMarker <- function(workflow_data, name) {
    if (!is.null(workflow_data[[name]])) return(workflow_data[[name]])
    if (protDiaWorkflowPayloadAvailable(workflow_data, name)) TRUE else NULL
}

resolveProtDiaWorkflowTable <- function(workflow_data, name) {
    if (!name %in% PROT_DIA_EVICT_FIELDS) {
        protDiaArtifactAbort(
            "DIA-NN payload resolver received an unsupported field",
            "multischolar_invalid_prot_dia_payload_role"
        )
    }
    resolveArtifactWorkflowTable(
        workflow_data,
        name,
        protDiaPayloadEvictionContract(),
        ref_fn = \(owner, field) {
            protDiaWorkflowPayloadRef(owner, field)
        },
        read_fn = protDiaArtifactReadTable
    )
}
