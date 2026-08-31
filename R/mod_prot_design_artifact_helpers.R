# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

protDiaArtifactMetadataTable <- function(value, owner) {
    artifactStageMetadataTable(value)
}

protDiaArtifactOptionalTable <- function(value, role) {
    artifactStageOptionalTable(value, role, protDiaArtifactAbort)
}

protDiaDesignArtifactTables <- function(workflow_data, state_object) {
    args <- if (isS4(state_object) && "args" %in% methods::slotNames(state_object)) {
        methods::slot(state_object, "args")
    } else {
        workflow_data$config_list
    }
    contrasts <- workflow_data$contrasts_tbl
    if (is.character(contrasts)) {
        contrasts <- data.frame(contrasts = contrasts, stringsAsFactors = FALSE)
    }
    list(
        cleaned_data = workflow_data$data_cln,
        design_matrix = workflow_data$design_matrix,
        contrasts = protDiaArtifactOptionalTable(contrasts, "contrasts"),
        args = protDiaArtifactMetadataTable(args, "args"),
        annotations = protDiaArtifactOptionalTable(
            workflow_data$uniprot_dat_cln,
            "annotations"
        ),
        sequences = protDiaArtifactOptionalTable(
            workflow_data$aa_seq_tbl_final,
            "sequences"
        )
    )
}

protDiaArtifactSetRunStatus <- function(
    connection,
    identity,
    stage,
    completed,
    timestamp
) {
    artifactStageSetRunStatus(
        connection,
        identity,
        stage,
        completed,
        timestamp,
        protDiaArtifactAbort
    )
}

protDiaArtifactCommitStageRefs <- function(
    session,
    connection,
    identity,
    stage,
    timestamp
) {
    artifactStageCommitRefs(
        session,
        connection,
        identity,
        stage,
        timestamp,
        protDiaArtifactAbort
    )
}

protDiaArtifactUpdateStageStatus <- function(
    context,
    stage,
    completed,
    failure_injector = NULL,
    resource_policy = NULL
) {
    artifactStageUpdateStatus(
        context,
        stage,
        completed,
        failure_injector,
        resource_policy,
        protDiaArtifactAbort
    )
}

protDiaArtifactStateAudit <- function(stage) {
    list(
        record_id = artifactOpaqueId("audit"),
        provenance_status = "artifact_dual_write_validated",
        stage_id = stage$stage_id,
        run_id = stage$run_id,
        action_id = stage$action_id,
        readthrough_contract_version = stage$readthrough_contract_version,
        parent_import_run_id = stage$parent_import$run_id,
        parent_import_generation_id = stage$parent_import$generation_id,
        parent_import_artifact_id = stage$parent_import$artifact_id,
        parent_import_semantic_digest = stage$parent_import$semantic_digest,
        stage_artifact_refs = stage$refs
    )
}

protDiaArtifactImportParent <- function(workflow_data) {
    result <- workflow_data$artifact_stage_results$import
    ref <- result$refs$canonical_data
    valid <- is.list(result) && isTRUE(result$enabled) && isTRUE(result$ok) &&
        isTRUE(result$committed) && workflowCapabilityScalarString(result$run_id) &&
        workflowCapabilityScalarString(result$generation_id) && is.list(ref) &&
        workflowCapabilityScalarString(ref$artifact_id) &&
        workflowCapabilityScalarString(ref$hash_policy$semantic$digest)
    if (!isTRUE(valid)) {
        return(list(
            run_id = NULL,
            generation_id = NULL,
            artifact_id = NULL,
            semantic_digest = NULL
        ))
    }
    list(
        run_id = result$run_id,
        generation_id = result$generation_id,
        artifact_id = ref$artifact_id,
        semantic_digest = ref$hash_policy$semantic$digest
    )
}

protDiaArtifactContrastsKind <- function(value) {
    artifactStageContrastsKind(value)
}

protDiaArtifactMemoryCheckpoint <- function(workflow_data, state_name) {
    manager <- workflow_data$state_manager
    state_object <- manager$getState(state_name)
    if (is.null(state_object) || !isS4(state_object)) {
        protDiaArtifactAbort(
            "DIA-NN design dual-write requires the completed memory S4 checkpoint",
            "multischolar_missing_prot_dia_memory_checkpoint"
        )
    }
    list(manager = manager, state_object = state_object)
}

protDiaArtifactWriteDesignStage <- function(
    prepared,
    workflow_data,
    state_name,
    checkpoint,
    failure_injector
) {
    parent_import <- protDiaArtifactImportParent(workflow_data)
    stage <- writeProtDiaStageArtifacts(
        prepared$context,
        stage_id = "design",
        tables = protDiaDesignArtifactTables(
            workflow_data,
            checkpoint$state_object
        ),
        parameters = list(
            state_name = state_name,
            workflow_type = workflowStateType(checkpoint$manager),
            readthrough_contract_version = 1L,
            parent_import_run_id = parent_import$run_id,
            parent_import_generation_id = parent_import$generation_id,
            parent_import_artifact_id = parent_import$artifact_id,
            parent_import_semantic_digest = parent_import$semantic_digest,
            contrasts_kind = protDiaArtifactContrastsKind(
                workflow_data$contrasts_tbl
            ),
            annotation_available = !is.null(workflow_data$uniprot_dat_cln),
            sequence_data_available = !is.null(workflow_data$aa_seq_tbl_final)
        ),
        deferred_commit = TRUE,
        failure_injector = failure_injector
    )
    stage$readthrough_contract_version <- 1L
    stage$parent_import <- parent_import
    stage
}

protDiaArtifactNewDesignStateManager <- function(prepared, manager_factory) {
    manager_factory(
        workflow_context = prepared$context,
        workflow_descriptor = artifactDiaWorkflowDescriptor(),
        descriptor_catalogue = artifactWorkflowDescriptorCatalogue(),
        codec_catalogue = artifactS4CodecCatalogue(),
        verify_hydration_fn = artifactWorkflowStateVerifyHydrationInline
    )
}

protDiaArtifactSaveDesignState <- function(
    manager,
    workflow_data,
    state_name,
    state_object,
    stage,
    save_state_fn
) {
    manager$setWorkflowType("DIA")
    audit <- protDiaArtifactStateAudit(stage)
    if (is.null(save_state_fn)) {
        manager$saveState(
            state_name = state_name,
            s4_data_object = state_object,
            config_object = workflow_data$config_list,
            description = "DIA-NN design dual-write checkpoint.",
            audit_metadata = audit
        )
    } else {
        save_state_fn(manager, state_name, state_object, audit)
    }
    hydration_verification <- manager$getLastHydrationVerification()
    deferred <- identical(
        hydration_verification$mode,
        "deferred_exact_pending"
    ) && is.null(hydration_verification$hydrated_digest)
    exact <- identical(
        hydration_verification$expected_digest,
        hydration_verification$hydrated_digest
    )
    if (!is.list(hydration_verification) ||
        !isTRUE(hydration_verification$valid) || (!deferred && !exact) ||
        isTRUE(hydration_verification$complete_payload_returned)) {
        protDiaArtifactAbort(
            "DIA-NN design artifact lacks an exact process-bound parity proof",
            "multischolar_inexact_prot_dia_s4_parity"
        )
    }
    list(
        state_manifest = manager$exportState(),
        state_metadata = manager$getStateMetadata(state_name),
        hydration_verification = hydration_verification
    )
}

protDiaArtifactFailDesignStage <- function(context, stage, state_error) {
    try(
        protDiaArtifactUpdateStageStatus(
            context,
            stage,
            completed = FALSE
        ),
        silent = TRUE
    )
    stop(state_error)
}

persistProtDiaDesignArtifactsCore <- function(
    workflow_data,
    state_name,
    failure_injector = NULL,
    manager_factory = newWorkflowState,
    save_state_fn = NULL
) {
    prepared <- prepareProtDiaArtifactContext(workflow_data)
    if (!isTRUE(prepared$enabled)) {
        return(list(
            enabled = FALSE,
            ok = TRUE,
            stage_id = "design",
            reason = prepared$reason
        ))
    }
    checkpoint <- protDiaArtifactMemoryCheckpoint(workflow_data, state_name)
    stage <- protDiaArtifactWriteDesignStage(
        prepared,
        workflow_data,
        state_name,
        checkpoint,
        failure_injector
    )
    manager <- NULL
    state_result <- tryCatch(
        {
            manager <- protDiaArtifactNewDesignStateManager(prepared, manager_factory)
            on.exit({
                if (!is.null(manager)) manager$close()
            }, add = TRUE)
            protDiaArtifactSaveDesignState(
                manager,
                workflow_data,
                state_name,
                checkpoint$state_object,
                stage,
                save_state_fn
            )
        },
        error = \(error) error
    )
    if (inherits(state_result, "error")) {
        if (!is.null(manager)) {
            try(manager$close(), silent = TRUE)
            manager <- NULL
        }
        protDiaArtifactFailDesignStage(prepared$context, stage, state_result)
    }
    manager$close()
    manager <- NULL
    protDiaArtifactUpdateStageStatus(
        prepared$context,
        stage,
        completed = TRUE,
        failure_injector = failure_injector
    )
    stage$committed <- TRUE
    c(list(enabled = TRUE, ok = TRUE), stage, state_result)
}

persistProtDiaDesignArtifacts <- function(
    workflow_data,
    state_name = workflowStateCurrentName(workflow_data$state_manager),
    failure_injector = NULL,
    manager_factory = newWorkflowState,
    save_state_fn = NULL,
    log_warn = logger::log_warn
) {
    runProtDiaArtifactSafely(
        workflow_data,
        "design",
        \() persistProtDiaDesignArtifactsCore(
            workflow_data,
            state_name,
            failure_injector = failure_injector,
            manager_factory = manager_factory,
            save_state_fn = save_state_fn
        ),
        log_warn = log_warn
    )
}
