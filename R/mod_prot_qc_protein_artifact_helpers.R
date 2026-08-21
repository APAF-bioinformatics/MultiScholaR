# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

.PROT_DIA_PROTEIN_QC_AUDIT_SCHEMA <- "multischolar.protein_qc_audit"
.PROT_DIA_PROTEIN_QC_AUDIT_VERSION <- 1L
.PROT_DIA_PROTEIN_QC_STAGES <- c(
    "protein_rollup",
    "protein_accession_cleanup",
    "protein_intensity_filter",
    "protein_duplicate_aggregation",
    "protein_replicate_filter"
)

protDiaProteinQcArtifactMode <- function() {
    mode <- getOption("multischolar.prot_dia.protein_qc_artifacts", "enabled")
    match.arg(mode, c("enabled", "disabled"))
}

protDiaProteinQcPackageVersion <- function(package) {
    tryCatch(
        as.character(utils::packageVersion(package)),
        error = \(error) "unavailable"
    )
}

protDiaProteinQcSoftware <- function(stage_id, parameters) {
    method <- parameters$rollup_method
    engine <- if (identical(method, "iq_maxlfq")) {
        list(name = "iq", version = protDiaProteinQcPackageVersion("iq"))
    } else if (identical(method, "limpa_dpc_quant")) {
        list(name = "limpa", version = protDiaProteinQcPackageVersion("limpa"))
    } else {
        list(name = "MultiScholaR", version = .peptideQcImplementationVersion())
    }
    list(
        name = "MultiScholaR",
        version = .peptideQcImplementationVersion(),
        source = "R package",
        stage_id = stage_id,
        scientific_engine = engine
    )
}

protDiaProteinQcTable <- function(object) {
    if (methods::is(object, "PeptideQuantitativeData")) {
        return(object@peptide_data)
    }
    if (methods::is(object, "ProteinQuantitativeData")) {
        return(object@protein_quant_table)
    }
    NULL
}

protDiaProteinQcActiveKey <- function(object) {
    if (methods::is(object, "PeptideQuantitativeData")) {
        return(object@protein_id_column)
    }
    if (methods::is(object, "ProteinQuantitativeData")) {
        return(resolveProteinQuantIdentityColumn(object))
    }
    NULL
}

protDiaProteinQcArgs <- function(object) {
    if (!isS4(object) || !"args" %in% methods::slotNames(object)) return(list())
    args <- object@args
    if (is.list(args)) args else list()
}

protDiaProteinQcCurrentRecord <- function(object) {
    audit <- protDiaProteinQcArgs(object)$protein_qc_audit
    if (!is.list(audit) || !is.list(audit$records) ||
        length(audit$records) == 0L) {
        return(NULL)
    }
    record <- tail(audit$records, 1L)[[1L]]
    if (!identical(record$record_id, audit$current_record_id)) {
        protDiaArtifactAbort(
            "DIA-NN protein-QC audit head does not match its latest record",
            "multischolar_invalid_prot_dia_protein_qc_audit"
        )
    }
    record
}

protDiaProteinQcPeptideAuditId <- function(object) {
    audit <- protDiaProteinQcArgs(object)$peptide_qc_audit
    record_id <- audit$current_record_id
    if (workflowStateScalarString(record_id)) record_id else "legacy_untracked_peptide"
}

protDiaProteinQcDigestArg <- function(object, name) {
    value <- protDiaProteinQcArgs(object)[[name]]
    if (is.null(value)) return("absent")
    .peptideQcDigest(value)
}

protDiaProteinQcIdTableDigest <- function(object) {
    if (!methods::is(object, "ProteinQuantitativeData")) return("not_protein")
    .peptideQcDigest(object@protein_id_table)
}

protDiaProteinQcAuditSummary <- function(object) {
    data <- protDiaProteinQcTable(object)
    key <- protDiaProteinQcActiveKey(object)
    if (!is.data.frame(data) || !workflowStateScalarString(key) ||
        !key %in% names(data)) {
        return(list(class = class(object)[1L], rows = NA_integer_))
    }
    list(
        class = class(object)[1L],
        rows = as.integer(nrow(data)),
        columns = as.integer(ncol(data)),
        active_protein_key = key,
        distinct_proteins = as.integer(length(unique(data[[key]]))),
        data_digest = .peptideQcDigest(data)
    )
}

protDiaProteinQcEmitAudit <- function(
    before,
    after,
    stage_id,
    parent_state,
    current_state,
    parent_generation_id,
    resolved_parameters,
    status,
    transformation_type,
    software,
    now
) {
    if (!methods::is(after, "ProteinQuantitativeData")) return(after)
    after_args <- protDiaProteinQcArgs(after)
    audit <- after_args$protein_qc_audit
    records <- audit$records
    if (!is.list(records)) records <- list()
    prior <- protDiaProteinQcCurrentRecord(before)
    parent_record_id <- if (is.null(prior)) {
        protDiaProteinQcPeptideAuditId(before)
    } else {
        prior$record_id
    }
    substantive <- list(
        schema = .PROT_DIA_PROTEIN_QC_AUDIT_SCHEMA,
        schema_version = .PROT_DIA_PROTEIN_QC_AUDIT_VERSION,
        stage_id = stage_id,
        status = status,
        transformation_type = transformation_type,
        parent_state = parent_state,
        current_state = current_state,
        parent_generation_id = parent_generation_id,
        parent_record_id = parent_record_id,
        peptide_audit_record_id = protDiaProteinQcPeptideAuditId(after),
        peptide_feature_key_map_digest = protDiaProteinQcDigestArg(
            after,
            "peptide_feature_key_map"
        ),
        protein_accession_provenance_digest = protDiaProteinQcDigestArg(
            after,
            "protein_accession_provenance"
        ),
        protein_id_table_digest = protDiaProteinQcIdTableDigest(after),
        resolved_parameters = .peptideQcCanonicalise(resolved_parameters),
        software = software,
        before_summary = protDiaProteinQcAuditSummary(before),
        after_summary = protDiaProteinQcAuditSummary(after)
    )
    record_id <- paste0(
        "protein-qc:",
        substr(.peptideQcDigest(substantive), 1L, 24L)
    )
    record <- c(
        substantive,
        list(
            record_id = record_id,
            timestamp_utc = format(now, "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC"),
            canonical_digest = .peptideQcDigest(substantive)
        )
    )
    records[[length(records) + 1L]] <- record
    after@args$protein_qc_audit <- list(
        schema = .PROT_DIA_PROTEIN_QC_AUDIT_SCHEMA,
        schema_version = .PROT_DIA_PROTEIN_QC_AUDIT_VERSION,
        records = records,
        current_record_id = record_id
    )
    after
}

protDiaProteinQcTransition <- function(before, after, stage_id) {
    if (!stage_id %in% .PROT_DIA_PROTEIN_QC_STAGES ||
        !methods::is(after, "ProteinQuantitativeData")) {
        return(FALSE)
    }
    if (identical(stage_id, "protein_rollup")) {
        return(methods::is(before, "PeptideQuantitativeData"))
    }
    methods::is(before, "ProteinQuantitativeData")
}

protDiaProteinQcArtifactEligible <- function(
    workflow_data,
    state_manager,
    before,
    after,
    stage_id
) {
    if (!protDiaProteinQcTransition(before, after, stage_id)) return(FALSE)
    if (inherits(state_manager, "ArtifactWorkflowState")) {
        return(identical(workflowStateType(state_manager), "DIA"))
    }
    !identical(protDiaProteinQcArtifactMode(), "disabled") &&
        protDiaPeptideQcWorkflowData(workflow_data) &&
        protDiaArtifactCoordinatorOwned(workflow_data)
}

protDiaProteinQcValidateParent <- function(manager, before) {
    parent <- manager$getState()
    if (!identical(parent, before)) {
        protDiaArtifactAbort(
            "DIA-NN protein-QC parent differs from the active scientific state",
            "multischolar_prot_dia_protein_qc_parent_mismatch"
        )
    }
    invisible(parent)
}

protDiaProteinQcRejectedReasons <- function(before, after, key_columns, stage_id) {
    parent_keys <- artifactWorkflowStateEntityKeys(
        before@protein_quant_table,
        key_columns
    )
    child_keys <- artifactWorkflowStateEntityKeys(
        after@protein_quant_table,
        key_columns
    )
    rejected <- parent_keys[!parent_keys %in% child_keys]
    if (length(rejected) == 0L) return(character())
    stats::setNames(
        rep(paste0("removed_by_", stage_id), length(rejected)),
        rejected
    )
}

protDiaProteinQcSelectionHint <- function(
    before,
    after,
    stage_id,
    state_name,
    parent_state,
    audit_parameters,
    transformation_type,
    software
) {
    if (!identical(transformation_type, "filter") ||
        !methods::is(before, "ProteinQuantitativeData")) {
        return(NULL)
    }
    key_columns <- protDiaProteinQcActiveKey(before)
    if (!workflowStateScalarString(key_columns) ||
        !identical(key_columns, protDiaProteinQcActiveKey(after))) {
        return(NULL)
    }
    record <- protDiaProteinQcCurrentRecord(after)
    if (is.null(record)) {
        protDiaArtifactAbort(
            "DIA-NN protein-QC artifact persistence requires its audit record",
            "multischolar_missing_prot_dia_protein_qc_audit"
        )
    }
    hint <- tryCatch(
        newArtifactRowSelectionHint(
            slot_name = "protein_quant_table",
            key_columns = key_columns,
            method = stage_id,
            normalized_parameters = .peptideQcCanonicalise(audit_parameters),
            software = software[c("name", "version", "source")],
            lineage = list(
                audit_record_id = record$record_id,
                state_name = state_name,
                parent_state = parent_state,
                parent_record_id = record$parent_record_id
            ),
            rejected_reasons = protDiaProteinQcRejectedReasons(
                before,
                after,
                key_columns,
                stage_id
            )
        ),
        multischolar_invalid_artifact_row_selection = \(error) NULL,
        multischolar_ambiguous_artifact_row_selection = \(error) NULL
    )
    if (is.null(hint)) return(NULL)
    tryCatch(
        {
            artifactWorkflowStateSelectionPlan(before, after, hint)
            hint
        },
        multischolar_artifact_selection_requires_materialization = \(error) NULL,
        multischolar_invalid_artifact_row_selection = \(error) NULL,
        multischolar_ambiguous_artifact_row_selection = \(error) NULL
    )
}

protDiaProteinQcSafeArtifactConfig <- function(config_object, software) {
    config <- config_object
    file_fields <- names(config)[grepl("_file$", names(config))]
    owned_files <- list()
    for (field in file_fields) {
        value <- config[[field]]
        if (workflowStateScalarString(value)) {
            owned_files[[field]] <- basename(value)
            config[[field]] <- basename(value)
        }
    }
    config$artifact_provenance <- list(
        software = software,
        output_ownership = list(
            scope = "workflow_session",
            files = owned_files
        )
    )
    config
}

protDiaProteinQcRecordResult <- function(workflow_data, result) {
    if (protDiaPeptideQcWorkflowData(workflow_data)) {
        recordProtDiaArtifactResult(workflow_data, result$stage_id, result)
    }
    invisible(result)
}

protDiaProteinQcSaveMemory <- function(
    state_manager,
    state_name,
    object,
    config_object,
    description,
    audit_metadata
) {
    save_state <- state_manager$saveState
    supported <- names(formals(save_state))
    args <- list(
        state_name = state_name,
        s4_data_object = object,
        config_object = config_object,
        description = description
    )
    if ("audit_metadata" %in% supported || "..." %in% supported) {
        args$audit_metadata <- audit_metadata
    }
    do.call(save_state, args)
}

saveProtDiaProteinQcArtifactState <- function(
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
    now,
    failure_injector = NULL
) {
    if (!protDiaProteinQcArtifactEligible(
        workflow_data,
        state_manager,
        before,
        after,
        stage_id
    )) {
        return(list(handled = FALSE, enabled = FALSE, object = after))
    }
    resources <- protDiaPeptideQcArtifactManager(workflow_data, state_manager)
    manager <- resources$manager
    if (isTRUE(resources$owned)) on.exit(manager$close(), add = TRUE)
    protDiaProteinQcValidateParent(manager, before)
    parent_state <- manager$getCurrentStateName()
    memory_parent_state <- workflowStateCurrentName(state_manager)
    if (isTRUE(resources$owned) && !identical(parent_state, memory_parent_state)) {
        protDiaArtifactAbort(
            "DIA-NN protein-QC managers disagree on the active parent state",
            "multischolar_prot_dia_protein_qc_parent_mismatch"
        )
    }
    parent_generation_id <- manager$getCurrentGenerationId()
    software <- protDiaProteinQcSoftware(stage_id, audit_parameters)
    after <- protDiaProteinQcEmitAudit(
        before,
        after,
        stage_id,
        parent_state,
        state_name,
        parent_generation_id,
        audit_parameters,
        status,
        transformation_type,
        software,
        now
    )
    record <- protDiaProteinQcCurrentRecord(after)
    config_object$audit_record_id <- record$record_id
    artifact_config <- protDiaProteinQcSafeArtifactConfig(config_object, software)
    hint <- protDiaProteinQcSelectionHint(
        before,
        after,
        stage_id,
        state_name,
        parent_state,
        audit_parameters,
        transformation_type,
        software
    )
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
        config_object = artifact_config,
        description = description,
        audit_metadata = record,
        persistence_hint = hint,
        failure_injector = failure_injector,
        action_id = artifactOpaqueId("action"),
        expected_parent_generation_id = parent_generation_id
    )
    if (!identical(result$status, "accepted")) {
        protDiaArtifactAbort(
            "DIA-NN protein-QC state commit did not advance its exact parent",
            "multischolar_prot_dia_protein_qc_commit_rejected",
            result = result
        )
    }
    if (isTRUE(resources$owned)) {
        memory_error <- tryCatch(
            {
                protDiaProteinQcSaveMemory(
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
                protDiaArtifactAbort(
                    paste(
                        "DIA-NN protein-QC compatibility save and artifact",
                        "rollback both failed"
                    ),
                    "multischolar_prot_dia_protein_qc_divergent_state",
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
            "DIA-NN protein-QC artifact hydration changed the scientific state",
            "multischolar_inexact_prot_dia_protein_qc_hydration"
        )
    }
    output <- list(
        handled = TRUE,
        enabled = TRUE,
        ok = TRUE,
        committed = TRUE,
        object = after,
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
    protDiaProteinQcRecordResult(workflow_data, output)
}

saveProtProteinQcState <- function(
    workflow_data,
    state_manager,
    before,
    after,
    stage_id,
    state_name,
    config_object,
    description,
    audit_parameters = config_object,
    status = "applied",
    transformation_type = "filter",
    now = Sys.time(),
    failure_injector = NULL
) {
    artifact_result <- saveProtDiaProteinQcArtifactState(
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
        now,
        failure_injector
    )
    if (isTRUE(artifact_result$handled)) return(artifact_result$object)
    parent_state <- workflowStateCurrentName(state_manager)
    software <- protDiaProteinQcSoftware(stage_id, audit_parameters)
    after <- protDiaProteinQcEmitAudit(
        before,
        after,
        stage_id,
        parent_state,
        state_name,
        paste0("memory:", parent_state),
        audit_parameters,
        status,
        transformation_type,
        software,
        now
    )
    record <- protDiaProteinQcCurrentRecord(after)
    config_object$audit_record_id <- record$record_id
    protDiaProteinQcSaveMemory(
        state_manager,
        state_name,
        after,
        config_object,
        description,
        record
    )
    after
}

protDiaProteinQcDualManagerMode <- function(workflow_data, state_manager) {
    !inherits(state_manager, "ArtifactWorkflowState") &&
        !identical(protDiaProteinQcArtifactMode(), "disabled") &&
        protDiaPeptideQcWorkflowData(workflow_data) &&
        protDiaArtifactCoordinatorOwned(workflow_data)
}

revertProtDiaProteinQcState <- function(workflow_data, state_name) {
    state_manager <- workflow_data$state_manager
    if (inherits(state_manager, "ArtifactWorkflowState")) {
        return(state_manager$revertToState(state_name))
    }
    if (!protDiaProteinQcDualManagerMode(workflow_data, state_manager)) {
        return(state_manager$revertToState(state_name))
    }
    artifact_manager <- newProtDiaArtifactStateManager(
        workflow_data$workflow_context
    )
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
        protDiaArtifactAbort(
            "DIA-NN protein-QC managers differ before revert",
            "multischolar_prot_dia_protein_qc_parent_mismatch"
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
        protDiaArtifactAbort(
            "DIA-NN protein-QC revert hydrated different active states",
            "multischolar_inexact_prot_dia_protein_qc_hydration"
        )
    }
    protDiaProteinQcRecordResult(workflow_data, list(
        enabled = TRUE,
        ok = TRUE,
        stage_id = "protein_qc_revert",
        state_name = state_name,
        artifact_generation_id = artifact_manager$getCurrentGenerationId()
    ))
    reverted_memory
}

initializeProtDiaProteinQcSessionContext <- function(
    workflow_data,
    experiment_paths,
    omic_type,
    experiment_label
) {
    if (!protDiaPeptideQcWorkflowData(workflow_data) ||
        !protDiaArtifactCoordinatorOwned(workflow_data)) {
        return(invisible(NULL))
    }
    existing <- workflow_data$protein_qc_context
    progress_env <- existing$progress_env
    if (!is.environment(progress_env)) progress_env <- new.env(parent = emptyenv())
    project_key <- paste0(omic_type, "_", experiment_label)
    workflow_data$protein_qc_context <- list(
        experiment_paths = experiment_paths,
        project_dirs = stats::setNames(list(experiment_paths), project_key),
        omic_type = omic_type,
        experiment_label = experiment_label,
        progress_env = progress_env
    )
    invisible(workflow_data$protein_qc_context)
}

protDiaProteinQcUpdateFiltering <- function(
    update_fn,
    workflow_data,
    data,
    step_name,
    omic_type,
    experiment_label,
    return_grid,
    overwrite
) {
    args <- list(
        data = data,
        step_name = step_name,
        omic_type = omic_type,
        experiment_label = experiment_label,
        return_grid = return_grid,
        overwrite = overwrite
    )
    context <- if (protDiaPeptideQcWorkflowData(workflow_data)) {
        workflow_data$protein_qc_context
    } else {
        NULL
    }
    if (protDiaArtifactCoordinatorOwned(workflow_data) && is.list(context) &&
        is.environment(context$progress_env)) {
        args$project_dirs <- context$project_dirs
        args$progress_env <- context$progress_env
    }
    do.call(update_fn, args)
}

protDiaProteinQcInvokeOutputs <- function(update_outputs_fn, args, workflow_data) {
    supported <- names(formals(update_outputs_fn))
    if ("workflowData" %in% supported || "..." %in% supported) {
        args$workflowData <- workflow_data
    }
    do.call(update_outputs_fn, args)
}
