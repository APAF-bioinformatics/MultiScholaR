# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

protDiaPeptideQcArtifactMode <- function() {
    mode <- getOption("multischolar.prot_dia.peptide_qc_artifacts", "enabled")
    match.arg(mode, c("enabled", "disabled"))
}

.PROT_DIA_PEPTIDE_QC_STAGES <- c(
    "qvalue_filter",
    "intensity_filter",
    "sample_filter",
    "protein_evidence_filter",
    "replicate_filter",
    "precursor_rollup",
    "imputation"
)

protDiaPeptideQcWorkflowData <- function(workflow_data) {
    is.environment(workflow_data) || inherits(workflow_data, "reactivevalues")
}

protDiaPeptideQcTransition <- function(before, after, stage_id) {
    methods::is(before, "PeptideQuantitativeData") &&
        methods::is(after, "PeptideQuantitativeData") &&
        stage_id %in% .PROT_DIA_PEPTIDE_QC_STAGES
}

protDiaPeptideQcStableKeyColumns <- function(object) {
    data <- .peptideQcObjectTable(object)
    roles <- .peptideQcColumnRoles(object)
    candidates <- unique(c(
        roles$active_protein_id,
        roles$run,
        roles$stripped_peptide,
        roles$modified_peptidoform,
        "Precursor.Id",
        "Precursor.Charge"
    ))
    candidates <- candidates[
        !is.na(candidates) & candidates %in% names(data)
    ]
    if (length(candidates) == 0L) return(NULL)
    valid <- tryCatch({
        artifactWorkflowStateEntityKeys(data, candidates)
        TRUE
    }, error = function(...) FALSE)
    if (isTRUE(valid)) candidates else NULL
}

protDiaPeptideQcCurrentAuditRecord <- function(object) {
    audit <- .peptideQcObjectArgs(object)$peptide_qc_audit
    if (!is.list(audit) || !is.list(audit$records) ||
        length(audit$records) == 0L) {
        return(NULL)
    }
    record <- tail(audit$records, 1L)[[1L]]
    if (!identical(record$record_id, audit$current_record_id)) {
        protDiaArtifactAbort(
            "DIA-NN peptide-QC audit head does not match its latest record",
            "multischolar_invalid_prot_dia_peptide_qc_audit"
        )
    }
    record
}

protDiaPeptideQcRejectedReasons <- function(before, after, key_columns) {
    parent_keys <- artifactWorkflowStateEntityKeys(
        .peptideQcObjectTable(before),
        key_columns
    )
    child_keys <- artifactWorkflowStateEntityKeys(
        .peptideQcObjectTable(after),
        key_columns
    )
    rejected_keys <- parent_keys[!parent_keys %in% child_keys]
    if (length(rejected_keys) == 0L) return(character())
    record <- protDiaPeptideQcCurrentAuditRecord(after)
    ledger <- record$removal_ledger
    if (!is.data.frame(ledger) || nrow(ledger) == 0L ||
        !all(key_columns %in% names(ledger))) {
        return(stats::setNames(
            rep("removed_by_stage_contract", length(rejected_keys)),
            rejected_keys
        ))
    }
    keys <- artifactWorkflowStateEntityKeys(ledger, key_columns)
    reasons <- if ("failure_reason" %in% names(ledger)) {
        as.character(ledger$failure_reason)
    } else {
        rep("removed_by_stage_contract", nrow(ledger))
    }
    ledger_reasons <- stats::setNames(reasons, keys)
    matched <- match(rejected_keys, names(ledger_reasons))
    rejected_reasons <- unname(ledger_reasons[matched])
    rejected_reasons[is.na(rejected_reasons) | !nzchar(rejected_reasons)] <-
        "removed_by_stage_contract"
    stats::setNames(rejected_reasons, rejected_keys)
}

protDiaPeptideQcSelectionHint <- function(
    before,
    after,
    stage_id,
    state_name,
    audit_parameters,
    status,
    transformation_type
) {
    selection_eligible <- methods::is(after, "PeptideQuantitativeData") &&
        (identical(transformation_type, "filter") ||
            status %in% c("skipped", "no_op"))
    if (!isTRUE(selection_eligible)) return(NULL)
    key_columns <- protDiaPeptideQcStableKeyColumns(before)
    if (is.null(key_columns) || !all(key_columns %in% names(after@peptide_data))) {
        return(NULL)
    }
    record <- protDiaPeptideQcCurrentAuditRecord(after)
    if (is.null(record)) {
        protDiaArtifactAbort(
            "DIA-NN peptide-QC artifact persistence requires its audit record",
            "multischolar_missing_prot_dia_peptide_qc_audit"
        )
    }
    hint <- newArtifactRowSelectionHint(
        slot_name = "peptide_data",
        key_columns = key_columns,
        method = stage_id,
        normalized_parameters = .peptideQcCanonicalise(audit_parameters),
        software = list(
            name = "MultiScholaR",
            version = .peptideQcImplementationVersion(),
            source = "R package"
        ),
        lineage = list(
            audit_record_id = record$record_id,
            state_name = state_name,
            parent_state = record$parent_state,
            parent_record_id = record$parent_record_id
        ),
        rejected_reasons = protDiaPeptideQcRejectedReasons(
            before,
            after,
            key_columns
        )
    )
    tryCatch({
        artifactWorkflowStateSelectionPlan(before, after, hint)
        hint
    },
    multischolar_artifact_selection_requires_materialization = \(error) NULL,
    multischolar_ambiguous_artifact_row_selection = \(error) NULL,
    multischolar_invalid_artifact_row_selection = \(error) NULL)
}

protDiaPeptideQcArtifactEligible <- function(
    workflow_data,
    state_manager,
    before,
    after,
    stage_id
) {
    if (!protDiaPeptideQcTransition(before, after, stage_id)) return(FALSE)
    if (inherits(state_manager, "ArtifactWorkflowState")) {
        return(identical(workflowStateType(state_manager), "DIA"))
    }
    if (identical(protDiaPeptideQcArtifactMode(), "disabled") ||
        !protDiaPeptideQcWorkflowData(workflow_data)) {
        return(FALSE)
    }
    protDiaArtifactCoordinatorOwned(workflow_data)
}

protDiaPeptideQcArtifactManager <- function(workflow_data, state_manager) {
    if (inherits(state_manager, "ArtifactWorkflowState")) {
        return(list(manager = state_manager, owned = FALSE))
    }
    list(
        manager = newProtDiaArtifactStateManager(
            workflow_data$workflow_context
        ),
        owned = TRUE
    )
}

protDiaPeptideQcValidateParent <- function(manager, before) {
    parent <- manager$getState()
    if (!isS4(parent) || !identical(class(parent), class(before))) {
        valid <- FALSE
    } else {
        scientific_slots <- setdiff(methods::slotNames(parent), "args")
        valid <- all(vapply(scientific_slots, \(slot_name) {
            identical(
                methods::slot(parent, slot_name),
                methods::slot(before, slot_name)
            )
        }, logical(1))) && identical(
            .peptideQcObjectArgs(parent)$peptide_qc_audit,
            .peptideQcObjectArgs(before)$peptide_qc_audit
        )
    }
    if (!isTRUE(valid)) {
        protDiaArtifactAbort(
            "DIA-NN peptide-QC parent differs from the active scientific state",
            "multischolar_prot_dia_peptide_qc_parent_mismatch"
        )
    }
    invisible(parent)
}

protDiaPeptideQcRecordArtifactResult <- function(workflow_data, result) {
    if (protDiaPeptideQcWorkflowData(workflow_data)) {
        recordProtDiaArtifactResult(workflow_data, result$stage_id, result)
    }
    invisible(result)
}

saveProtDiaPeptideQcArtifactState <- function(
    workflow_data,
    state_manager,
    before,
    after,
    stage_id,
    state_name,
    config_object,
    description,
    audit_parameters,
    status,
    transformation_type,
    failure_injector = NULL
) {
    if (!protDiaPeptideQcArtifactEligible(
        workflow_data,
        state_manager,
        before,
        after,
        stage_id
    )) {
        return(list(handled = FALSE, enabled = FALSE))
    }
    resources <- protDiaPeptideQcArtifactManager(workflow_data, state_manager)
    manager <- resources$manager
    if (isTRUE(resources$owned)) {
        on.exit(manager$close(), add = TRUE)
    }
    protDiaPeptideQcValidateParent(manager, before)
    parent_state <- manager$getCurrentStateName()
    memory_parent_state <- .peptideQcCurrentStateName(state_manager)
    if (isTRUE(resources$owned) && !identical(parent_state, memory_parent_state)) {
        protDiaArtifactAbort(
            "DIA-NN peptide-QC managers disagree on the active parent state",
            "multischolar_prot_dia_peptide_qc_parent_mismatch"
        )
    }
    parent_generation_id <- manager$getCurrentGenerationId()
    hint <- protDiaPeptideQcSelectionHint(
        before,
        after,
        stage_id,
        state_name,
        audit_parameters,
        status,
        transformation_type
    )
    record <- protDiaPeptideQcCurrentAuditRecord(after)
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
        config_object = config_object,
        description = description,
        audit_metadata = record,
        persistence_hint = hint,
        failure_injector = failure_injector,
        action_id = artifactOpaqueId("action"),
        expected_parent_generation_id = parent_generation_id
    )
    if (!identical(result$status, "accepted")) {
        protDiaArtifactAbort(
            "DIA-NN peptide-QC state commit did not advance its exact parent",
            "multischolar_prot_dia_peptide_qc_commit_rejected",
            result = result
        )
    }
    if (isTRUE(resources$owned)) {
        memory_error <- tryCatch({
            state_manager$saveState(
                state_name = state_name,
                s4_data_object = after,
                config_object = config_object,
                description = description
            )
            NULL
        }, error = identity)
        if (inherits(memory_error, "error")) {
            restore_error <- tryCatch({
                manager$revertToState(parent_state)
                NULL
            }, error = identity)
            if (inherits(restore_error, "error")) {
                protDiaArtifactAbort(
                    paste(
                        "DIA-NN peptide-QC compatibility save and artifact",
                        "rollback both failed"
                    ),
                    "multischolar_prot_dia_peptide_qc_divergent_state",
                    memory_error = memory_error,
                    restore_error = restore_error
                )
            }
            stop(memory_error)
        }
    }
    hydrated <- manager$getState(state_name)
    if (!identical(hydrated, after)) {
        protDiaArtifactAbort(
            "DIA-NN peptide-QC artifact hydration changed the scientific state",
            "multischolar_inexact_prot_dia_peptide_qc_hydration"
        )
    }
    output <- list(
        handled = TRUE,
        enabled = TRUE,
        ok = TRUE,
        committed = TRUE,
        stage_id = stage_id,
        state_name = state_name,
        generation_id = result$generation_id,
        parent_generation_id = parent_generation_id,
        representation = result$representation,
        metrics = result$metrics,
        audit_record_id = record$record_id,
        compaction = list(
            enabled = FALSE,
            reason = "representative_measurement_not_available"
        )
    )
    protDiaPeptideQcRecordArtifactResult(workflow_data, output)
}

protDiaPeptideQcDualManagerMode <- function(workflow_data, state_manager) {
    !inherits(state_manager, "ArtifactWorkflowState") &&
        !identical(protDiaPeptideQcArtifactMode(), "disabled") &&
        protDiaPeptideQcWorkflowData(workflow_data) &&
        protDiaArtifactCoordinatorOwned(workflow_data)
}

protDiaPeptideQcRestoreRevert <- function(
    artifact_manager,
    state_manager,
    artifact_state,
    memory_state
) {
    errors <- list()
    errors$artifact <- tryCatch({
        artifact_manager$revertToState(artifact_state)
        NULL
    }, error = identity)
    errors$memory <- tryCatch({
        state_manager$revertToState(memory_state)
        NULL
    }, error = identity)
    errors <- Filter(\(error) inherits(error, "error"), errors)
    if (length(errors) > 0L) {
        protDiaArtifactAbort(
            "DIA-NN peptide-QC revert rollback could not restore both managers",
            "multischolar_prot_dia_peptide_qc_divergent_state",
            restore_errors = errors
        )
    }
    invisible(TRUE)
}

revertProtDiaPeptideQcState <- function(workflow_data, state_name) {
    state_manager <- workflow_data$state_manager
    if (inherits(state_manager, "ArtifactWorkflowState")) {
        return(state_manager$revertToState(state_name))
    }
    if (!protDiaPeptideQcDualManagerMode(workflow_data, state_manager)) {
        return(state_manager$revertToState(state_name))
    }

    artifact_manager <- newProtDiaArtifactStateManager(
        workflow_data$workflow_context
    )
    on.exit(artifact_manager$close(), add = TRUE)
    memory_current <- state_manager$getCurrentStateName()
    artifact_current <- artifact_manager$getCurrentStateName()
    memory_object <- state_manager$getState()
    artifact_object <- artifact_manager$getState()
    targets_exist <- isTRUE(state_manager$hasState(state_name)) &&
        isTRUE(artifact_manager$hasState(state_name))
    if (!targets_exist || !identical(memory_current, artifact_current) ||
        !identical(memory_object, artifact_object)) {
        protDiaArtifactAbort(
            "DIA-NN peptide-QC managers differ before revert",
            "multischolar_prot_dia_peptide_qc_parent_mismatch"
        )
    }
    memory_target <- state_manager$getState(state_name)
    artifact_target <- artifact_manager$getState(state_name)
    if (!identical(memory_target, artifact_target)) {
        protDiaArtifactAbort(
            "DIA-NN peptide-QC revert targets differ between managers",
            "multischolar_prot_dia_peptide_qc_parent_mismatch"
        )
    }

    reverted_artifact <- artifact_manager$revertToState(state_name)
    memory_error <- tryCatch({
        state_manager$revertToState(state_name)
        NULL
    }, error = identity)
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
        protDiaArtifactAbort(
            "DIA-NN peptide-QC revert hydrated different active states",
            "multischolar_inexact_prot_dia_peptide_qc_hydration"
        )
    }
    protDiaPeptideQcRecordArtifactResult(workflow_data, list(
        enabled = TRUE,
        ok = TRUE,
        stage_id = "peptide_qc_revert",
        state_name = state_name,
        artifact_generation_id = artifact_manager$getCurrentGenerationId()
    ))
    reverted_memory
}
