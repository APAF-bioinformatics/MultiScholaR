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
    context <- workflow_data$workflow_context
    if (!inherits(context, "WorkflowContext") || !context$isBound()) return(NULL)
    rollout_fn(context)
}

protDiaEvictionProofFlags <- function(proof) {
    names <- c(
        "payload_validated", "registry_ready", "current_pointer_verified",
        "readthrough_completed", "annotation_completed"
    )
    stats::setNames(
        vapply(names, \(name) isTRUE(proof[[name]]), logical(1)),
        names
    )
}

protDiaEvictionCheckpointReady <- function(checkpoint, proof) {
    is.list(checkpoint) && isTRUE(checkpoint$verified) &&
        identical(checkpoint$strategy, "read_only_artifact_reconstruction") &&
        identical(checkpoint$project_id, proof$project_id) &&
        identical(checkpoint$workflow_id, proof$workflow_id) &&
        identical(checkpoint$design_run_id, proof$design_run_id) &&
        identical(
            checkpoint$state_generation_id,
            proof$state_generation_id
        ) && identical(
            checkpoint$consumer_inventory_digest,
            protDiaEvictionInventoryDigest()
        )
}

protDiaEvictionCacheReady <- function(workflow_state) {
    if (!inherits(workflow_state, "ArtifactWorkflowState")) return(TRUE)
    info <- workflow_state$getCacheInfo()
    is.list(info) && is.numeric(info$entries) && info$entries <= 1L
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
    proof <- workflow_data$artifact_readthrough_proof
    checkpoint <- workflow_data$artifact_compatibility_checkpoint
    source_ready <- vapply(
        PROT_DIA_EVICT_FIELDS,
        \(name) !is.null(workflow_data[[name]]),
        logical(1)
    )
    c(
        exact_dia_canary = protDiaArtifactCoordinatorOwned(workflow_data),
        evict_rollout = identical(
            protDiaEvictionRollout(workflow_data, rollout_fn),
            "evict"
        ),
        source_payload_available = all(source_ready),
        protDiaEvictionProofFlags(proof),
        compatibility_checkpoint_verified = protDiaEvictionCheckpointReady(
            checkpoint,
            proof
        ),
        consumer_inventory_verified = identical(
            proof$consumer_inventory_digest,
            protDiaEvictionInventoryDigest()
        ),
        single_cache_entry = protDiaEvictionCacheReady(
            workflow_data$state_manager
        )
    )
}

validateProtDiaEvictionReadiness <- function(readiness) {
    required <- c(
        "exact_dia_canary", "evict_rollout", "source_payload_available",
        "payload_validated", "registry_ready", "current_pointer_verified",
        "readthrough_completed", "annotation_completed",
        "compatibility_checkpoint_verified", "consumer_inventory_verified",
        "single_cache_entry"
    )
    valid <- is.logical(readiness) && identical(names(readiness), required) &&
        !anyNA(readiness)
    if (!isTRUE(valid)) {
        protDiaArtifactAbort(
            "DIA-NN eviction readiness is malformed",
            "multischolar_invalid_prot_dia_eviction_readiness"
        )
    }
    readiness
}

protDiaEvictionSourceFieldBytes <- function(workflow_data) {
    vapply(
        PROT_DIA_EVICT_FIELDS,
        \(name) as.numeric(utils::object.size(workflow_data[[name]])),
        numeric(1)
    )
}

setProtDiaWorkflowField <- function(workflow_data, name, value) {
    workflow_data[[name]] <- value
    invisible(TRUE)
}

clearProtDiaWorkflowPayloads <- function(
    workflow_data,
    clear_fn,
    release_cache_fn
) {
    original <- lapply(PROT_DIA_EVICT_FIELDS, \(name) workflow_data[[name]])
    names(original) <- PROT_DIA_EVICT_FIELDS
    committed <- FALSE
    on.exit({
        if (!committed) {
            for (name in PROT_DIA_EVICT_FIELDS) {
                setProtDiaWorkflowField(workflow_data, name, original[[name]])
            }
        }
    }, add = TRUE)
    for (name in PROT_DIA_EVICT_FIELDS) {
        clear_fn(workflow_data, name, NULL)
    }
    release_cache_fn(workflow_data$state_manager)
    cleared <- vapply(
        PROT_DIA_EVICT_FIELDS,
        \(name) is.null(workflow_data[[name]]),
        logical(1)
    )
    if (!all(cleared)) {
        protDiaArtifactAbort(
            "DIA-NN reactive payload eviction did not clear every source",
            "multischolar_incomplete_prot_dia_eviction"
        )
    }
    committed <- TRUE
    invisible(TRUE)
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
    gc_fn = \() gc(full = TRUE)
) {
    previous <- workflow_data$artifact_stage_results$eviction
    if (is.list(previous) && isTRUE(previous$evicted)) return(previous)
    readiness <- validateProtDiaEvictionReadiness(
        readiness_fn(workflow_data, rollout_fn)
    )
    failed <- names(readiness)[!readiness]
    if (length(failed) > 0L) {
        return(list(
            enabled = isTRUE(readiness[["exact_dia_canary"]]) &&
                isTRUE(readiness[["evict_rollout"]]),
            ok = TRUE,
            evicted = FALSE,
            reason = "eviction_prerequisites_incomplete",
            failed_prerequisites = failed
        ))
    }
    source_field_bytes <- protDiaEvictionSourceFieldBytes(workflow_data)
    clearProtDiaWorkflowPayloads(
        workflow_data,
        clear_fn,
        release_cache_fn
    )
    invisible(gc_fn())
    result <- list(
        enabled = TRUE,
        ok = TRUE,
        evicted = TRUE,
        reason = "artifact_payloads_evicted",
        evicted_fields = PROT_DIA_EVICT_FIELDS,
        released_source_field_bytes = source_field_bytes,
        released_source_bytes_upper_bound = sum(source_field_bytes),
        state_generation_id = workflow_data$artifact_readthrough_proof$state_generation_id,
        compatibility_strategy = "read_only_artifact_reconstruction"
    )
    recordProtDiaArtifactResult(workflow_data, "eviction", result)
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
    resource_policy = NULL
) {
    prepared <- createProtDiaResumeContext(
        experiment_paths,
        experiment_label,
        storage_policy
    )
    if (!isTRUE(prepared$enabled)) return(prepared)
    bundle <- hydrateProtDiaResumeBundle(
        prepared$context,
        resource_policy,
        retain_source_payloads = FALSE
    )
    on.exit(bundle$state_manager$close(), add = TRUE)
    workflow_data$artifact_readthrough_refs <- list(
        import = bundle$evidence$import$refs,
        design = bundle$evidence$design$refs
    )
    readthrough <- recordProtDiaReadthroughProof(workflow_data, bundle)
    list(
        enabled = TRUE,
        ok = TRUE,
        resumed = FALSE,
        verified = TRUE,
        reason = "artifact_readthrough_verified_for_eviction",
        project_id = bundle$evidence$identity$project_id,
        import_run_id = bundle$evidence$import$run_id,
        design_run_id = bundle$evidence$design$run_id,
        state_generation_id = bundle$evidence$current_state$generation_id,
        compatibility_checkpoint = readthrough$compatibility_checkpoint
    )
}

settleProtDiaArtifactWorkflowSafely <- function(
    workflow_data,
    experiment_paths,
    experiment_label,
    storage_policy = NULL,
    resource_policy = NULL,
    rollout_fn = \(context) context$getStorageDecision()$effective_rollout,
    log_warn = logger::log_warn
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
                    resource_policy
                )
            } else {
                verification <- list(verified = TRUE)
            }
            if (!isTRUE(verification$verified)) {
                verification
            } else {
                evictProtDiaWorkflowPayloads(
                    workflow_data,
                    rollout_fn = rollout_fn
                )
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

protDiaWorkflowPayloadAvailable <- function(workflow_data, name) {
    if (!name %in% PROT_DIA_EVICT_FIELDS) return(!is.null(workflow_data[[name]]))
    if (!is.null(workflow_data[[name]])) return(TRUE)
    refs <- workflow_data$artifact_readthrough_refs
    if (identical(name, "data_tbl")) {
        return(is.list(refs$import$canonical_data))
    }
    is.list(refs$design$cleaned_data)
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
    if (!is.null(workflow_data[[name]])) return(workflow_data[[name]])
    if (!protDiaArtifactCoordinatorOwned(workflow_data)) return(NULL)
    refs <- workflow_data$artifact_readthrough_refs
    ref <- if (identical(name, "data_tbl")) {
        refs$import$canonical_data
    } else {
        refs$design$cleaned_data
    }
    if (!is.list(ref)) return(NULL)
    identity <- workflow_data$workflow_context$getIdentity()
    store <- newArtifactStore(
        workflow_data$workflow_context$getPaths(),
        identity$project_id
    )
    protDiaArtifactReadTable(store, ref)
}
