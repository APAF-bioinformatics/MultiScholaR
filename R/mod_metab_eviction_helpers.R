# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

METAB_EVICT_FIELDS <- c("data_tbl", "data_cln")

metabEvictionStageGate <- function() {
    list(
        minimum_released_source_bytes = 1,
        minimum_source_graph_reduction = 0.25,
        maximum_post_eviction_cache_entries = 0L
    )
}

metabEvictionConsumerInventory <- function() {
    inventory <- metabReadthroughConsumerInventory()
    inventory$descriptor_id <-
        artifactMetabolomicsWorkflowDescriptor()$descriptor_id
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

metabPayloadEvictionContract <- function() {
    newArtifactPayloadEvictionContract(
        descriptor = artifactMetabolomicsWorkflowDescriptor(),
        owner_label = "metabolomics",
        owner_condition_name = "exact_metabolomics_tuple",
        source_fields = METAB_EVICT_FIELDS,
        inventory_digest_fn = \() artifactSemanticDigest(
            metabReadthroughConsumerInventory()
        ),
        compatibility_strategy =
            "restore_memory_reads_keep_metabolomics_generations",
        abort_fn = metabArtifactAbort,
        invalid_readiness_class =
            "multischolar_invalid_metabolomics_eviction_readiness",
        incomplete_eviction_class =
            "multischolar_incomplete_metabolomics_eviction"
    )
}

metabEvictionReadiness <- function(
    workflow_data,
    rollout_fn = \(context) context$getStorageDecision()$effective_rollout
) {
    artifactPayloadEvictionReadiness(
        workflow_data,
        metabPayloadEvictionContract(),
        rollout_fn
    )
}

evictMetabArtifactWorkflowPayloads <- function(
    workflow_data,
    rollout_fn = \(context) context$getStorageDecision()$effective_rollout,
    readiness_fn = NULL,
    clear_fn = setArtifactWorkflowField,
    release_cache_fn = workflowStateReleaseHydrationCache,
    gc_fn = \() gc(full = TRUE)
) {
    result <- evictArtifactWorkflowPayloads(
        workflow_data,
        metabPayloadEvictionContract(),
        rollout_fn,
        readiness_fn,
        clear_fn,
        release_cache_fn
    )
    if (isTRUE(result$evicted)) invisible(gc_fn())
    result
}

evictMetabArtifactWorkflowPayloadsSafely <- function(
    workflow_data,
    rollout_fn = \(context) context$getStorageDecision()$effective_rollout,
    log_warn = logger::log_warn,
    ...
) {
    result <- tryCatch(
        evictMetabArtifactWorkflowPayloads(
            workflow_data,
            rollout_fn,
            ...
        ),
        error = \(error) {
            log_warn(paste(
                "metabolomics eviction failed; source values were restored:",
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

metabWorkflowPayloadRefs <- function(workflow_data, name) {
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

resolveMetabWorkflowAssays <- function(workflow_data, name) {
    if (!name %in% METAB_EVICT_FIELDS) {
        metabArtifactAbort(
            "metabolomics payload resolver received an unsupported field",
            "multischolar_invalid_metabolomics_eviction_readiness"
        )
    }
    if (!is.null(workflow_data[[name]])) return(workflow_data[[name]])
    if (!metabArtifactCoordinatorOwned(workflow_data)) return(NULL)
    refs <- metabWorkflowPayloadRefs(workflow_data, name)
    if (!is.list(refs) || length(refs) == 0L) return(NULL)
    context <- workflow_data$workflow_context
    store <- newArtifactStore(
        context$getPaths(),
        context$getIdentity()$project_id
    )
    values <- lapply(refs, \(ref) artifactStageReadTable(
        metabArtifactReadthroughAdapter(),
        store,
        ref
    ))
    names(values) <- names(refs)
    values
}

metabWorkflowPayloadAvailable <- function(workflow_data, name) {
    if (!name %in% METAB_EVICT_FIELDS) return(!is.null(workflow_data[[name]]))
    !is.null(workflow_data[[name]]) ||
        length(metabWorkflowPayloadRefs(workflow_data, name)) > 0L
}

evaluateMetabEvictionStageGate <- function(
    before_source_graph_bytes,
    eviction_result,
    state_manager
) {
    gate <- metabEvictionStageGate()
    reduction <- if (before_source_graph_bytes > 0) 1 else 0
    cache_entries <- if (inherits(state_manager, "ArtifactWorkflowState")) {
        state_manager$getCacheInfo()$entries
    } else {
        0L
    }
    checks <- c(
        source_bytes = isTRUE(
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

settleMetabArtifactWorkflowSafely <- function(
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
    evictMetabArtifactWorkflowPayloadsSafely(
        workflow_data,
        rollout_fn,
        log_warn
    )
}
