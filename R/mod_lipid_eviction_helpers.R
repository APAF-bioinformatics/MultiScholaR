# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

LIPID_EVICT_FIELDS <- c("data_tbl", "data_cln")

lipidEvictionStageGate <- function() {
    list(
        minimum_released_source_bytes = 1,
        minimum_source_graph_reduction = 0.25,
        maximum_post_eviction_cache_entries = 0L
    )
}

lipidEvictionConsumerInventory <- function() {
    inventory <- lipidReadthroughConsumerInventory()
    inventory$descriptor_id <-
        artifactLipidomicsWorkflowDescriptor()$descriptor_id
    inventory$source_after_eviction <- c(
        "ordered import refs", "design refs and exact S4", "exact S4 state",
        "exact S4 state", "exact S4 plus coordinator metadata",
        "exact S4/design/contrasts", "exact S4 plus coordinator metadata"
    )
    inventory[c(
        "descriptor_id", "category", "reader", "restored_source",
        "source_after_eviction", "verified"
    )]
}

lipidPayloadEvictionContract <- function() {
    newArtifactPayloadEvictionContract(
        descriptor = artifactLipidomicsWorkflowDescriptor(),
        owner_label = "lipidomics",
        owner_condition_name = "exact_lipidomics_tuple",
        source_fields = LIPID_EVICT_FIELDS,
        inventory_digest_fn = \() artifactSemanticDigest(
            lipidReadthroughConsumerInventory()
        ),
        compatibility_strategy =
            "restore_memory_reads_keep_lipidomics_generations",
        abort_fn = lipidArtifactAbort,
        invalid_readiness_class =
            "multischolar_invalid_lipidomics_eviction_readiness",
        incomplete_eviction_class =
            "multischolar_incomplete_lipidomics_eviction"
    )
}

lipidEvictionReadiness <- function(
    workflow_data,
    rollout_fn = \(context) context$getStorageDecision()$effective_rollout
) {
    artifactPayloadEvictionReadiness(
        workflow_data,
        lipidPayloadEvictionContract(),
        rollout_fn
    )
}

evictLipidArtifactWorkflowPayloads <- function(
    workflow_data,
    rollout_fn = \(context) context$getStorageDecision()$effective_rollout,
    readiness_fn = NULL,
    clear_fn = setArtifactWorkflowField,
    release_cache_fn = workflowStateReleaseHydrationCache,
    gc_fn = \() gc(full = TRUE)
) {
    result <- evictArtifactWorkflowPayloads(
        workflow_data,
        lipidPayloadEvictionContract(),
        rollout_fn,
        readiness_fn,
        clear_fn,
        release_cache_fn
    )
    if (isTRUE(result$evicted)) {
        result$registry_session_released <- isTRUE(
            workflowStateReleaseRegistrySession(workflow_data$state_manager)
        )
        invisible(gc_fn())
        result$allocator_pages_released <- isTRUE(
            artifactReleaseProcessAllocator()
        )
    }
    result
}

evictLipidArtifactWorkflowPayloadsSafely <- function(
    workflow_data,
    rollout_fn = \(context) context$getStorageDecision()$effective_rollout,
    log_warn = logger::log_warn,
    ...
) {
    result <- tryCatch(
        evictLipidArtifactWorkflowPayloads(
            workflow_data,
            rollout_fn,
            ...
        ),
        error = \(error) {
            log_warn(paste(
                "lipidomics eviction failed; source values were restored:",
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

lipidWorkflowPayloadRefs <- function(workflow_data, name) {
    proof <- workflow_data$artifact_readthrough_proof
    assay_order <- proof$assay_order
    refs <- workflow_data$artifact_readthrough_refs
    if (!is.character(assay_order) || length(assay_order) == 0L) return(NULL)
    if (name == "data_tbl") {
        roles <- sprintf("assay_%04d", seq_along(assay_order))
        selected <- refs$import[roles]
    } else {
        roles <- sprintf("cleaned_assay_%04d", seq_along(assay_order))
        selected <- refs$design[roles]
    }
    names(selected) <- assay_order
    selected
}

resolveLipidWorkflowAssays <- function(workflow_data, name) {
    if (!name %in% LIPID_EVICT_FIELDS) {
        lipidArtifactAbort(
            "lipidomics payload resolver received an unsupported field",
            "multischolar_invalid_lipidomics_eviction_readiness"
        )
    }
    if (!is.null(workflow_data[[name]])) return(workflow_data[[name]])
    if (!lipidArtifactCoordinatorOwned(workflow_data)) return(NULL)
    refs <- lipidWorkflowPayloadRefs(workflow_data, name)
    if (!is.list(refs) || length(refs) == 0L) return(NULL)
    context <- workflow_data$workflow_context
    store <- newArtifactStore(
        context$getPaths(),
        context$getIdentity()$project_id
    )
    values <- lapply(refs, \(ref) artifactSettledResumeReadLocator(
        lipidArtifactReadthroughAdapter(),
        store,
        ref
    ))
    names(values) <- names(refs)
    values
}

lipidWorkflowPayloadAvailable <- function(workflow_data, name) {
    if (!name %in% LIPID_EVICT_FIELDS) return(!is.null(workflow_data[[name]]))
    !is.null(workflow_data[[name]]) ||
        length(lipidWorkflowPayloadRefs(workflow_data, name)) > 0L
}

lipidWorkflowAssayNames <- function(workflow_data, name = "data_tbl") {
    value <- workflow_data[[name]]
    if (!is.null(value)) return(names(value))
    refs <- lipidWorkflowPayloadRefs(workflow_data, name)
    if (is.null(refs)) character() else names(refs)
}

evaluateLipidEvictionStageGate <- function(
    before_source_graph_bytes,
    eviction_result,
    state_manager
) {
    gate <- lipidEvictionStageGate()
    reduction <- if (before_source_graph_bytes > 0) 1 else 0
    cache_entries <- if (inherits(state_manager, "ArtifactWorkflowState")) {
        state_manager$getCacheInfo()$entries
    } else {
        0L
    }
    checks <- c(
        source_bytes = isTRUE(eviction_result$source_hydration_avoided) || isTRUE(
            eviction_result$released_source_bytes_upper_bound >=
                gate$minimum_released_source_bytes
        ),
        source_graph = reduction >= gate$minimum_source_graph_reduction,
        cache = cache_entries <= gate$maximum_post_eviction_cache_entries
    )
    list(
        passed = all(checks),
        checks = checks,
        measured = list(
            before_source_graph_bytes = before_source_graph_bytes,
            after_source_graph_bytes = 0,
            source_graph_reduction = reduction,
            released_source_bytes =
                eviction_result$released_source_bytes_upper_bound,
            post_eviction_cache_entries = cache_entries
        ),
        gate = gate
    )
}

settleLipidArtifactWorkflowSafely <- function(
    workflow_data,
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
    settled <- settledArtifactWorkflowPayloadResult(
        workflow_data,
        lipidPayloadEvictionContract()
    )
    if (!is.null(settled)) return(settled)
    evictLipidArtifactWorkflowPayloadsSafely(
        workflow_data,
        rollout_fn,
        log_warn
    )
}

settlePersistedLipidArtifactWorkflowSafely <- function(
    workflow_data,
    log_warn = logger::log_warn
) {
    rollout <- artifactPayloadEvictionRollout(
        workflow_data,
        \(context) context$getStorageDecision()$effective_rollout
    )
    if (!identical(rollout, "evict")) {
        return(list(
            enabled = FALSE,
            ok = TRUE,
            evicted = FALSE,
            reason = "rollout_below_evict"
        ))
    }
    tryCatch({
        bundle <- hydrateLipidResumeBundle(
            workflow_data$workflow_context,
            retain_source_payloads = FALSE
        )
        applyLipidResumeBundle(workflow_data, bundle)
        settleLipidArtifactWorkflowSafely(workflow_data, log_warn = log_warn)
    }, error = \(error) {
        log_warn(paste(
            "lipidomics post-design settlement failed; memory state retained:",
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
    })
}
