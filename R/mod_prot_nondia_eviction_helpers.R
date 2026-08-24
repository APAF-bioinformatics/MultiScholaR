# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

PROT_NONDIA_EVICT_FIELDS <- c("data_tbl", "data_cln")

#' Define independent non-DIA implementation-stage eviction gates
#' @return A named list of exact tuple stage gates.
#' @noRd
protNonDiaEvictionStageGates <- function() {
    gate <- list(
        minimum_released_source_bytes = 1,
        minimum_source_graph_reduction = 0.25,
        maximum_post_eviction_cache_entries = 0L
    )
    ids <- names(artifactProteomicsNonDiaDescriptorSpecs())
    stats::setNames(rep(list(gate), length(ids)), ids)
}

#' Resolve one tuple-specific non-DIA eviction stage gate
#' @param capability_id Exact workflow capability identifier.
#' @return The exact tuple implementation-stage gate.
#' @noRd
protNonDiaEvictionStageGate <- function(capability_id) {
    gate <- protNonDiaEvictionStageGates()[[capability_id]]
    if (is.null(gate)) {
        protNonDiaArtifactAbort(
            "non-DIA eviction has no exact tuple stage gate",
            "multischolar_missing_prot_nondia_eviction_gate",
            capability_id = capability_id
        )
    }
    gate
}

#' Inventory non-DIA eviction consumers for one exact tuple
#' @param descriptor Exact supported workflow descriptor.
#' @return A tuple-bound consumer inventory data frame.
#' @noRd
protNonDiaEvictionConsumerInventory <- function(descriptor) {
    descriptor <- protNonDiaReadthroughDescriptor(descriptor$descriptor_id)
    inventory <- protNonDiaReadthroughConsumerInventory()
    inventory$descriptor_id <- descriptor$descriptor_id
    inventory$workflow_type <- descriptor$identity$workflow_type
    inventory$input_format <- descriptor$identity$input_format
    inventory[, c(
        "descriptor_id", "workflow_type", "input_format", "category",
        "reader", "restored_source", "verified"
    )]
}

#' Digest one tuple-specific non-DIA eviction consumer inventory
#' @param descriptor Exact supported workflow descriptor.
#' @return A semantic digest.
#' @noRd
protNonDiaEvictionInventoryDigest <- function(descriptor) {
    artifactSemanticDigest(protNonDiaEvictionConsumerInventory(descriptor))
}

#' Define one exact non-DIA payload eviction contract
#' @param descriptor Exact supported workflow descriptor.
#' @return A descriptor-bound payload eviction contract.
#' @noRd
protNonDiaPayloadEvictionContract <- function(descriptor) {
    descriptor <- protNonDiaReadthroughDescriptor(descriptor$descriptor_id)
    newArtifactPayloadEvictionContract(
        descriptor = descriptor,
        owner_label = "non-DIA proteomics",
        owner_condition_name = "exact_nondia_tuple",
        source_fields = PROT_NONDIA_EVICT_FIELDS,
        inventory_digest_fn = \() {
            proof_inventory <- protNonDiaReadthroughConsumerInventory()
            artifactSemanticDigest(proof_inventory)
        },
        compatibility_strategy =
            "restore_memory_reads_keep_artifact_generations",
        abort_fn = protNonDiaArtifactAbort,
        invalid_readiness_class =
            "multischolar_invalid_prot_nondia_eviction_readiness",
        incomplete_eviction_class =
            "multischolar_incomplete_prot_nondia_eviction"
    )
}

#' Evaluate exact non-DIA payload eviction readiness
#' @param workflow_data Mutable proteomics workflow state.
#' @param descriptor Exact supported workflow descriptor.
#' @param rollout_fn Effective rollout provider.
#' @return A named logical readiness vector.
#' @noRd
protNonDiaEvictionReadiness <- function(
    workflow_data,
    descriptor,
    rollout_fn = \(context) context$getStorageDecision()$effective_rollout
) {
    artifactPayloadEvictionReadiness(
        workflow_data,
        protNonDiaPayloadEvictionContract(descriptor),
        rollout_fn
    )
}

#' Evict exact non-DIA import/design source payloads
#' @param workflow_data Mutable proteomics workflow state.
#' @param descriptor Exact supported workflow descriptor.
#' @param rollout_fn Effective rollout provider.
#' @param readiness_fn Optional readiness provider.
#' @param clear_fn Injectable workflow field writer.
#' @param release_cache_fn Injectable hydration-cache release function.
#' @return A payload-free eviction result.
#' @noRd
evictProtNonDiaWorkflowPayloads <- function(
    workflow_data,
    descriptor,
    rollout_fn = \(context) context$getStorageDecision()$effective_rollout,
    readiness_fn = NULL,
    clear_fn = setArtifactWorkflowField,
    release_cache_fn = workflowStateReleaseHydrationCache
) {
    if (is.null(readiness_fn)) {
        readiness_fn <- \(owner, resolver) {
            protNonDiaEvictionReadiness(owner, descriptor, resolver)
        }
    }
    evictArtifactWorkflowPayloads(
        workflow_data,
        protNonDiaPayloadEvictionContract(descriptor),
        rollout_fn,
        readiness_fn,
        clear_fn,
        release_cache_fn
    )
}

#' Safely evict exact non-DIA import/design source payloads
#' @param workflow_data Mutable proteomics workflow state.
#' @param descriptor Exact supported workflow descriptor.
#' @param rollout_fn Effective rollout provider.
#' @param log_warn Warning logger used for additive eviction failures.
#' @param ... Additional eviction controls used by tests.
#' @return The recorded eviction result, invisibly.
#' @noRd
evictProtNonDiaWorkflowPayloadsSafely <- function(
    workflow_data,
    descriptor,
    rollout_fn = \(context) context$getStorageDecision()$effective_rollout,
    log_warn = logger::log_warn,
    ...
) {
    result <- tryCatch(
        evictProtNonDiaWorkflowPayloads(
            workflow_data,
            descriptor,
            rollout_fn,
            ...
        ),
        error = \(error) {
            log_warn(paste(
                "non-DIA payload eviction failed; source values were restored:",
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
    recordArtifactStageResult(workflow_data, "eviction", result)
}

#' Verify settled non-DIA read-through before current-session eviction
#' @param workflow_data Mutable proteomics workflow state.
#' @param experiment_paths Workflow project paths.
#' @param experiment_label Project experiment label.
#' @param descriptor Exact supported workflow descriptor.
#' @param storage_policy Optional workflow storage policy.
#' @param resource_policy Optional project registry resource policy.
#' @return A payload-free settled read-through verification.
#' @noRd
verifyProtNonDiaReadthroughForEviction <- function(
    workflow_data,
    experiment_paths,
    experiment_label,
    descriptor,
    storage_policy = NULL,
    resource_policy = NULL
) {
    descriptor <- protNonDiaReadthroughDescriptor(descriptor$descriptor_id)
    prepared <- createProtNonDiaResumeContext(
        experiment_paths,
        experiment_label,
        descriptor$descriptor_id,
        storage_policy
    )
    if (!isTRUE(prepared$enabled)) return(prepared)
    bundle <- hydrateProtNonDiaResumeBundle(
        prepared$context,
        descriptor,
        resource_policy,
        retain_source_payloads = FALSE
    )
    on.exit(bundle$state_manager$close(), add = TRUE)
    workflow_data$artifact_readthrough_refs <- list(
        import = bundle$evidence$import$refs,
        design = bundle$evidence$design$refs
    )
    proof <- newProtNonDiaReadthroughProof(bundle)
    proof$consumer_inventory_digest <- artifactSemanticDigest(
        protNonDiaReadthroughConsumerInventory()
    )
    checkpoint <- newProtNonDiaCompatibilityCheckpoint(proof)
    workflow_data$artifact_readthrough_proof <- proof
    workflow_data$artifact_compatibility_checkpoint <- checkpoint
    list(
        enabled = TRUE,
        ok = TRUE,
        resumed = FALSE,
        verified = TRUE,
        reason = "artifact_readthrough_verified_for_eviction",
        capability_id = descriptor$descriptor_id,
        project_id = bundle$evidence$identity$project_id,
        import_run_id = bundle$evidence$import$run_id,
        design_run_id = bundle$evidence$design$run_id,
        state_generation_id = bundle$evidence$current_state$generation_id,
        compatibility_checkpoint = checkpoint
    )
}

#' Settle one current non-DIA workflow when eviction is certified
#' @param workflow_data Mutable proteomics workflow state.
#' @param experiment_paths Workflow project paths.
#' @param experiment_label Project experiment label.
#' @param descriptor Exact supported workflow descriptor.
#' @param storage_policy Optional workflow storage policy.
#' @param resource_policy Optional project registry resource policy.
#' @param rollout_fn Effective rollout provider.
#' @param log_warn Warning logger used for additive settlement failures.
#' @return The recorded settlement result, invisibly.
#' @noRd
settleProtNonDiaArtifactWorkflowSafely <- function(
    workflow_data,
    experiment_paths,
    experiment_label,
    descriptor,
    storage_policy = NULL,
    resource_policy = NULL,
    rollout_fn = \(context) context$getStorageDecision()$effective_rollout,
    log_warn = logger::log_warn
) {
    rollout <- artifactPayloadEvictionRollout(workflow_data, rollout_fn)
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
                verification <- verifyProtNonDiaReadthroughForEviction(
                    workflow_data,
                    experiment_paths,
                    experiment_label,
                    descriptor,
                    storage_policy,
                    resource_policy
                )
            } else {
                verification <- list(verified = TRUE)
            }
            if (!isTRUE(verification$verified)) {
                verification
            } else {
                evictProtNonDiaWorkflowPayloads(
                    workflow_data,
                    descriptor,
                    rollout_fn
                )
            }
        },
        error = \(error) {
            log_warn(paste(
                "non-DIA artifact settlement failed without evicting sources:",
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
    recordArtifactStageResult(workflow_data, "eviction", result)
}

#' Test whether a non-DIA workflow payload is available
#' @param workflow_data Mutable proteomics workflow state.
#' @param name Workflow field name.
#' @return A scalar logical.
#' @noRd
protNonDiaWorkflowPayloadAvailable <- function(workflow_data, name) {
    if (!name %in% PROT_NONDIA_EVICT_FIELDS) {
        return(!is.null(workflow_data[[name]]))
    }
    if (!is.null(workflow_data[[name]])) return(TRUE)
    refs <- workflow_data$artifact_readthrough_refs
    if (identical(name, "data_tbl")) {
        return(is.list(refs$import$canonical_data))
    }
    is.list(refs$design$cleaned_data)
}

#' Resolve one explicitly requested non-DIA source table
#' @param workflow_data Mutable proteomics workflow state.
#' @param name Exact source field name.
#' @return The retained or independently reconstructed table.
#' @noRd
resolveProtNonDiaWorkflowTable <- function(workflow_data, name) {
    if (!inherits(workflow_data$workflow_context, "WorkflowContext") ||
        !workflow_data$workflow_context$isBound()) {
        return(NULL)
    }
    descriptor <- findArtifactWorkflowDescriptor(
        workflow_data$workflow_context$getIdentity(),
        artifactWorkflowDescriptorCatalogue()
    )
    descriptor <- protNonDiaReadthroughDescriptor(descriptor$descriptor_id)
    adapter <- protNonDiaArtifactReadthroughAdapter(descriptor)
    resolveArtifactWorkflowTable(
        workflow_data,
        name,
        protNonDiaPayloadEvictionContract(descriptor),
        ref_fn = \(owner, field) {
            refs <- owner$artifact_readthrough_refs
            if (identical(field, "data_tbl")) {
                refs$import$canonical_data
            } else {
                refs$design$cleaned_data
            }
        },
        read_fn = \(store, ref) artifactStageReadTable(adapter, store, ref)
    )
}

#' Evaluate one tuple's implementation-stage memory gate
#' @param capability_id Exact workflow capability identifier.
#' @param before_source_graph_bytes Reachable source graph before eviction.
#' @param eviction_result Successful eviction result.
#' @param state_manager Current artifact state manager.
#' @return A structured gate decision and measured values.
#' @noRd
evaluateProtNonDiaEvictionStageGate <- function(
    capability_id,
    before_source_graph_bytes,
    eviction_result,
    state_manager
) {
    gate <- protNonDiaEvictionStageGate(capability_id)
    after_source_graph_bytes <- 0
    reduction <- if (before_source_graph_bytes > 0) {
        1 - after_source_graph_bytes / before_source_graph_bytes
    } else {
        0
    }
    cache_entries <- if (inherits(state_manager, "ArtifactWorkflowState")) {
        state_manager$getCacheInfo()$entries
    } else {
        0L
    }
    checks <- c(
        source_bytes = isTRUE(eviction_result$released_source_bytes_upper_bound >=
            gate$minimum_released_source_bytes),
        source_graph = reduction >= gate$minimum_source_graph_reduction,
        cache = cache_entries <= gate$maximum_post_eviction_cache_entries
    )
    list(
        capability_id = capability_id,
        passed = all(checks),
        checks = checks,
        measured = list(
            before_source_graph_bytes = before_source_graph_bytes,
            after_source_graph_bytes = after_source_graph_bytes,
            source_graph_reduction = reduction,
            released_source_bytes =
                eviction_result$released_source_bytes_upper_bound,
            post_eviction_cache_entries = cache_entries
        ),
        gate = gate
    )
}
