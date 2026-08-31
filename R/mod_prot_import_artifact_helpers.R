# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

protDiaArtifactAbort <- function(message, class, ...) {
    rlang::abort(
        message,
        class = c(class, "multischolar_prot_dia_artifact_error"),
        ...
    )
}

protDiaArtifactIdentityMatches <- function(identity) {
    artifactStageIdentityMatches(identity, artifactDiaWorkflowDescriptor()) ||
        artifactStageIdentityMatches(identity, artifactDiaWorkflowDescriptorV1())
}

protDiaArtifactCoordinatorOwned <- function(workflow_data) {
    context <- tryCatch(workflow_data$workflow_context, error = function(...) NULL)
    if (!inherits(context, "WorkflowContext") || !context$isBound()) return(FALSE)
    decision <- context$getStorageDecision()
    descriptor <- if (identical(decision$capability_version, "1.0.0")) {
        artifactDiaWorkflowDescriptorV1()
    } else {
        artifactDiaWorkflowDescriptor()
    }
    artifactStageCoordinatorOwned(workflow_data, descriptor)
}

prepareProtDiaArtifactContext <- function(
    workflow_data,
    format = workflow_data$data_format,
    data_type = workflow_data$data_type,
    descriptor_catalogue = artifactWorkflowDescriptorCatalogue()
) {
    if (!identical(format, "diann") || !identical(data_type, "peptide")) {
        return(list(enabled = FALSE, reason = "not_exact_dia_canary"))
    }
    prepareArtifactStageContext(
        workflow_data,
        workflow_type = "DIA",
        input_format = format,
        data_level = data_type,
        descriptor_catalogue = descriptor_catalogue
    )
}

recordProtDiaArtifactResult <- function(workflow_data, stage_id, result) {
    recordArtifactStageResult(workflow_data, stage_id, result)
}

runProtDiaArtifactSafely <- function(
    workflow_data,
    stage_id,
    operation,
    log_warn = logger::log_warn
) {
    runArtifactStageSafely(
        workflow_data,
        stage_id,
        operation,
        "DIA-NN",
        log_warn
    )
}

prepareProtDiaArtifactContextSafely <- function(
    workflow_data,
    format = workflow_data$data_format,
    data_type = workflow_data$data_type,
    log_warn = logger::log_warn
) {
    runProtDiaArtifactSafely(
        workflow_data,
        "context",
        \() {
            result <- prepareProtDiaArtifactContext(
                workflow_data,
                format = format,
                data_type = data_type
            )
            result$ok <- TRUE
            result$stage_id <- "context"
            result
        },
        log_warn = log_warn
    )
}

protDiaArtifactDescriptorContract <- function() {
    artifactStageDescriptorContract(artifactDiaWorkflowDescriptor())
}

protDiaArtifactPortableTable <- function(value, owner) {
    artifactStagePortableTable(value, owner, protDiaArtifactAbort)
}

protDiaArtifactPackageVersion <- function() {
    artifactStagePackageVersion()
}

protDiaArtifactSourceMetadata <- function(source_path, retain_source_uri = FALSE) {
    if (!artifactResourceScalarString(source_path) || !file.exists(source_path) ||
        dir.exists(source_path)) {
        protDiaArtifactAbort(
            "DIA-NN artifact source must be one existing file",
            "multischolar_invalid_prot_dia_source"
        )
    }
    extension <- tolower(tools::file_ext(source_path))
    list(
        source_role = "search_results",
        source_uri = if (isTRUE(retain_source_uri)) basename(source_path) else NULL,
        source_digest = artifactByteDigest(source_path),
        source_size_bytes = unname(as.numeric(file.info(source_path)$size)),
        parser_id = "MultiScholaR::importDIANNData",
        parser_version = "1.0.0",
        format_id = if (identical(extension, "parquet")) {
            "diann.parquet"
        } else {
            "diann.tsv"
        },
        data_level = "peptide"
    )
}

protDiaArtifactParameterJson <- function(value) {
    artifactStageParameterJson(value)
}

protDiaArtifactVerifyRef <- function(store, ref, expected) {
    managed <- artifactStoreManagedPaths(store, ref$logical_key, ref$artifact_id)
    sidecar <- artifactStoreReadSidecar(
        store,
        managed$sidecar,
        validate_payload = TRUE
    )
    payload <- arrow::read_parquet(
        artifactStoreResolveFile(store, ref$relative_path, must_exist = TRUE),
        as_data_frame = FALSE
    )
    hydrated <- decodeArtifactRectangular(payload, sidecar$codec_metadata)
    if (!identical(hydrated, expected)) {
        protDiaArtifactAbort(
            sprintf(
                "independent DIA-NN artifact hydration changed '%s'",
                ref$logical_key$state_role
            ),
            "multischolar_inexact_prot_dia_artifact_hydration"
        )
    }
    invisible(hydrated)
}

protDiaArtifactRegistryRun <- function(
    session,
    identity,
    run_id,
    action_id,
    status,
    timestamp
) {
    artifactStageRegistryRun(
        session,
        identity,
        run_id,
        action_id,
        status,
        timestamp
    )
}

protDiaArtifactRegistrySource <- function(
    session,
    identity,
    run_id,
    source,
    timestamp
) {
    artifactStageRegistrySource(session, identity, run_id, source, timestamp)
}

protDiaArtifactRegistryParameters <- function(
    session,
    identity,
    run_id,
    parameters,
    timestamp
) {
    artifactStageRegistryParameters(
        session,
        identity,
        run_id,
        parameters,
        timestamp
    )
}

protDiaArtifactRegistrySoftware <- function(session, identity, run_id, timestamp) {
    artifactStageRegistrySoftware(session, identity, run_id, timestamp)
}

protDiaArtifactRegistryRefs <- function(
    session,
    identity,
    run_id,
    refs,
    status,
    timestamp,
    register_run_artifacts
) {
    artifactStageRegistryRefs(
        session,
        identity,
        run_id,
        refs,
        status,
        timestamp,
        register_run_artifacts
    )
}

protDiaArtifactStageResources <- function(
    context,
    stage_id,
    resource_policy = NULL
) {
    if (!inherits(context, "WorkflowContext") || !context$isBound() ||
        !protDiaArtifactIdentityMatches(context$getIdentity())) {
        protDiaArtifactAbort(
            "DIA-NN stage persistence requires its exact bound context",
            "multischolar_invalid_prot_dia_artifact_context"
        )
    }
    artifactStageResources(
        context,
        artifactDiaWorkflowDescriptor(),
        stage_id,
        resource_policy,
        protDiaArtifactAbort
    )
}

protDiaArtifactWriteStageRefs <- function(
    stage,
    tables,
    failure_injector = NULL,
    verification = c("inline_exact", "deferred_exact")
) {
    artifactStageWriteRefs(
        stage,
        tables,
        failure_injector,
        verification,
        protDiaArtifactAbort
    )
}

protDiaArtifactRegisterStage <- function(
    stage,
    refs,
    parameters,
    source,
    deferred_commit
) {
    artifactStageRegister(stage, refs, parameters, source, deferred_commit)
}

writeProtDiaStageArtifacts <- function(
    context,
    stage_id,
    tables,
    parameters,
    source = NULL,
    deferred_commit = FALSE,
    failure_injector = NULL,
    resource_policy = NULL
) {
    writeArtifactStage(
        context,
        artifactDiaWorkflowDescriptor(),
        stage_id,
        tables,
        parameters,
        source,
        deferred_commit,
        failure_injector,
        resource_policy,
        protDiaArtifactAbort
    )
}

persistProtDiaImportArtifacts <- function(
    workflow_data,
    data_import_result,
    source_path,
    use_precursor_norm = TRUE,
    sanitize_names = FALSE,
    retain_source_uri = FALSE,
    failure_injector = NULL,
    pending_stage = NULL,
    worker_attempted = !is.null(pending_stage),
    log_warn = logger::log_warn
) {
    runProtDiaArtifactSafely(
        workflow_data,
        "import",
        \() {
            if (isTRUE(worker_attempted)) {
                if (is.null(pending_stage)) {
                    return(list(
                        enabled = TRUE,
                        ok = FALSE,
                        stage_id = "import",
                        reason = "artifact_worker_failed_no_retry",
                        committed = FALSE
                    ))
                }
                return(publishProtDiaArtifactWorkerStage(
                    workflow_data,
                    pending_stage
                ))
            }
            prepared <- prepareProtDiaArtifactContext(workflow_data)
            if (!isTRUE(prepared$enabled)) {
                return(list(
                    enabled = FALSE,
                    ok = TRUE,
                    stage_id = "import",
                    reason = prepared$reason
                ))
            }
            stage <- writeProtDiaStageArtifacts(
                prepared$context,
                stage_id = "import",
                tables = list(canonical_data = workflow_data$data_tbl),
                parameters = list(
                    column_mapping = data_import_result$column_mapping,
                    column_mapping_serialized = artifactWorkflowStateSerializeMetadata(
                        data_import_result$column_mapping,
                        "DIA-NN import column mapping"
                    ),
                    readthrough_contract_version = 1L,
                    use_precursor_norm = isTRUE(use_precursor_norm),
                    sanitize_names = isTRUE(sanitize_names)
                ),
                source = protDiaArtifactSourceMetadata(
                    source_path,
                    retain_source_uri = retain_source_uri
                ),
                failure_injector = failure_injector
            )
            c(list(enabled = TRUE, ok = TRUE), stage)
        },
        log_warn = log_warn
    )
}
