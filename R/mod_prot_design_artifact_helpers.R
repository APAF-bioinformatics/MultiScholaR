# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

protDiaArtifactMetadataTable <- function(value, owner) {
    if (is.null(value)) {
        return(data.frame(
            key = character(),
            value_json = character(),
            value_digest = character(),
            stringsAsFactors = FALSE
        ))
    }
    if (!is.list(value)) value <- list(value = value)
    value_names <- names(value)
    if (is.null(value_names) || any(!nzchar(value_names))) {
        value_names <- sprintf("item_%04d", seq_along(value))
    }
    data.frame(
        key = value_names,
        value_json = vapply(value, protDiaArtifactParameterJson, character(1)),
        value_digest = vapply(value, artifactSemanticDigest, character(1)),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
}

protDiaArtifactOptionalTable <- function(value, role) {
    if (is.null(value)) {
        return(data.frame(
            artifact_role = character(),
            availability = character(),
            stringsAsFactors = FALSE
        ))
    }
    protDiaArtifactPortableTable(value, role)
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
    run_status <- if (isTRUE(completed)) "completed" else "failed"
    run_count <- projectRegistryExecuteBound(
        connection,
        paste0(
            "UPDATE workflow_runs SET status = ?, completed_at = ?, updated_at = ? ",
            "WHERE project_id = ? AND workflow_id = ? AND run_id = ? ",
            "AND status = 'running'"
        ),
        list(
            run_status,
            timestamp,
            timestamp,
            identity$project_id,
            identity$workflow_id,
            stage$run_id
        )
    )
    if (!identical(as.integer(run_count), 1L)) {
        protDiaArtifactAbort(
            "DIA-NN artifact stage run status did not advance exactly once",
            "multischolar_prot_dia_registry_compare_and_set_failed"
        )
    }
    invisible(TRUE)
}

protDiaArtifactCommitStageRefs <- function(
    session,
    connection,
    identity,
    stage,
    timestamp
) {
    artifact_count <- projectRegistryExecuteBound(
        connection,
        paste0(
            "UPDATE artifacts SET status = 'committed', updated_at = ? ",
            "WHERE project_id = ? AND workflow_id = ? AND run_id = ? ",
            "AND status = 'validated'"
        ),
        list(
            timestamp,
            identity$project_id,
            identity$workflow_id,
            stage$run_id
        )
    )
    if (!identical(as.integer(artifact_count), as.integer(length(stage$refs)))) {
        protDiaArtifactAbort(
            "DIA-NN artifact stage refs did not commit as one complete set",
            "multischolar_prot_dia_registry_compare_and_set_failed"
        )
    }
    for (index in seq_along(stage$refs)) {
        ref <- stage$refs[[index]]
        projectRegistryWrite(session, "run_artifact", list(
            workflow_id = identity$workflow_id,
            run_id = stage$run_id,
            artifact_id = ref$artifact_id,
            direction = "output",
            artifact_role = ref$logical_key$state_role,
            ordinal = index - 1L,
            recorded_at = timestamp
        ))
    }
    invisible(TRUE)
}

protDiaArtifactUpdateStageStatus <- function(
    context,
    stage,
    completed,
    failure_injector = NULL,
    resource_policy = NULL
) {
    identity <- context$getIdentity()
    registry <- projectRegistryForContext(context, resource_policy)
    session <- initializeProjectRegistry(registry)
    on.exit(closeProjectRegistry(session), add = TRUE)
    connection <- projectRegistrySessionConnection(session, write = TRUE)
    timestamp <- artifactRefUtcNow()
    artifactStoreInvokeFailure(
        failure_injector,
        "before_stage_status_commit",
        list(stage_id = stage$stage_id, run_id = stage$run_id)
    )
    artifactWorkflowStateTransaction(session, \() {
        protDiaArtifactSetRunStatus(
            connection,
            identity,
            stage,
            completed,
            timestamp
        )
        if (isTRUE(completed)) {
            protDiaArtifactCommitStageRefs(
                session,
                connection,
                identity,
                stage,
                timestamp
            )
        }
    })
    invisible(TRUE)
}

protDiaArtifactStateAudit <- function(stage) {
    list(
        record_id = artifactOpaqueId("audit"),
        provenance_status = "artifact_dual_write_validated",
        stage_id = stage$stage_id,
        run_id = stage$run_id,
        action_id = stage$action_id,
        stage_artifact_refs = stage$refs
    )
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
    writeProtDiaStageArtifacts(
        prepared$context,
        stage_id = "design",
        tables = protDiaDesignArtifactTables(
            workflow_data,
            checkpoint$state_object
        ),
        parameters = list(
            state_name = state_name,
            workflow_type = workflowStateType(checkpoint$manager),
            annotation_available = !is.null(workflow_data$uniprot_dat_cln),
            sequence_data_available = !is.null(workflow_data$aa_seq_tbl_final)
        ),
        deferred_commit = TRUE,
        failure_injector = failure_injector
    )
}

protDiaArtifactNewDesignStateManager <- function(prepared, manager_factory) {
    manager_factory(
        workflow_context = prepared$context,
        workflow_descriptor = artifactDiaWorkflowDescriptor(),
        descriptor_catalogue = artifactWorkflowDescriptorCatalogue(),
        codec_catalogue = artifactS4CodecCatalogue()
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
    hydrated <- manager$getState(state_name)
    if (!identical(hydrated, state_object)) {
        protDiaArtifactAbort(
            "DIA-NN design artifact differs from its memory checkpoint",
            "multischolar_inexact_prot_dia_artifact_hydration"
        )
    }
    list(
        state_manifest = manager$exportState(),
        state_metadata = manager$getStateMetadata(state_name)
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
