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
    if (is.null(identity)) return(FALSE)
    expected <- artifactDiaWorkflowDescriptor()$identity
    identical(workflowCapabilityKey(identity), workflowCapabilityKey(expected))
}

protDiaArtifactCoordinatorOwned <- function(workflow_data) {
    context <- workflow_data$workflow_context
    if (!inherits(context, "WorkflowContext")) return(FALSE)
    if (!context$isBound()) return(FALSE)
    decision <- context$getStorageDecision()
    identical(decision$effective_backend, "artifact") &&
        protDiaArtifactIdentityMatches(context$getIdentity())
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
    context <- workflow_data$workflow_context
    if (!inherits(context, "WorkflowContext")) {
        return(list(enabled = FALSE, reason = "workflow_context_absent"))
    }
    if (!context$isBound()) {
        bindWorkflowContextFromImport(
            context,
            workflow_type = "DIA",
            input_format = "diann",
            data_level = "peptide",
            descriptor_catalogue = descriptor_catalogue
        )
    }
    decision <- context$getStorageDecision()
    enabled <- identical(decision$effective_backend, "artifact") &&
        identical(decision$effective_rollout, "dual_write") &&
        protDiaArtifactIdentityMatches(context$getIdentity())
    list(
        enabled = enabled,
        reason = decision$reason_code,
        context = context,
        decision = unclass(decision)
    )
}

recordProtDiaArtifactResult <- function(workflow_data, stage_id, result) {
    results <- workflow_data$artifact_stage_results
    if (is.null(results)) results <- list()
    results[[stage_id]] <- result
    workflow_data$artifact_stage_results <- results
    invisible(result)
}

runProtDiaArtifactSafely <- function(
    workflow_data,
    stage_id,
    operation,
    log_warn = logger::log_warn
) {
    result <- tryCatch(
        operation(),
        error = \(error) {
            log_warn(paste(
                "DIA-NN artifact",
                stage_id,
                "dual-write failed without changing memory state:",
                conditionMessage(error)
            ))
            list(
                enabled = TRUE,
                ok = FALSE,
                stage_id = stage_id,
                error_class = class(error),
                error_message = conditionMessage(error)
            )
        }
    )
    recordProtDiaArtifactResult(workflow_data, stage_id, result)
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
    descriptor <- artifactDiaWorkflowDescriptor()
    list(
        descriptor_id = descriptor$descriptor_id,
        descriptor_version = descriptor$descriptor_version,
        descriptor_digest = descriptor$descriptor_digest
    )
}

protDiaArtifactPortableTable <- function(value, owner) {
    if (!is.data.frame(value)) {
        protDiaArtifactAbort(
            sprintf("DIA-NN artifact '%s' must be a data frame", owner),
            "multischolar_invalid_prot_dia_artifact_table"
        )
    }
    if (inherits(value, "spec_tbl_df")) value <- tibble::as_tibble(value)
    value
}

protDiaArtifactPackageVersion <- function() {
    version <- tryCatch(
        as.character(utils::packageVersion("MultiScholaR")),
        error = \(...) NA_character_
    )
    if (!is.na(version)) return(version)
    description <- tryCatch(read.dcf("DESCRIPTION"), error = \(...) NULL)
    if (is.null(description)) "unknown" else description[[1L, "Version"]]
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
    as.character(jsonlite::toJSON(
        value,
        auto_unbox = TRUE,
        null = "null",
        na = "string",
        digits = NA
    ))
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
    projectRegistryWrite(session, "run", list(
        workflow_id = identity$workflow_id,
        run_id = run_id,
        status = status,
        action_id = action_id,
        started_at = timestamp,
        completed_at = if (identical(status, "completed")) timestamp else NULL,
        created_at = timestamp,
        updated_at = timestamp
    ))
}

protDiaArtifactRegistrySource <- function(
    session,
    identity,
    run_id,
    source,
    timestamp
) {
    if (is.null(source)) return(invisible(FALSE))
    projectRegistryWrite(session, "source", c(list(
        workflow_id = identity$workflow_id,
        run_id = run_id,
        source_id = artifactOpaqueId("source")
    ), source, list(recorded_at = timestamp)))
}

protDiaArtifactRegistryParameters <- function(
    session,
    identity,
    run_id,
    parameters,
    timestamp
) {
    for (parameter_name in names(parameters)) {
        value <- parameters[[parameter_name]]
        projectRegistryWrite(session, "parameter", list(
            workflow_id = identity$workflow_id,
            run_id = run_id,
            parameter_id = artifactOpaqueId("parameter"),
            parameter_name = parameter_name,
            value_json = protDiaArtifactParameterJson(value),
            value_digest = artifactSemanticDigest(value),
            recorded_at = timestamp
        ))
    }
    invisible(TRUE)
}

protDiaArtifactRegistrySoftware <- function(session, identity, run_id, timestamp) {
    version <- protDiaArtifactPackageVersion()
    projectRegistryWrite(session, "software", list(
        workflow_id = identity$workflow_id,
        run_id = run_id,
        software_id = artifactOpaqueId("software"),
        software_name = "MultiScholaR",
        software_version = version,
        software_source = "R package",
        software_digest = artifactSemanticDigest(list(
            package = "MultiScholaR",
            version = version
        )),
        recorded_at = timestamp
    ))
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
    for (index in seq_along(refs)) {
        record <- artifactWorkflowStateArtifactRecord(
            identity,
            refs[[index]],
            index - 1L
        )
        record$run_id <- run_id
        record$status <- status
        record$updated_at <- timestamp
        projectRegistryWrite(session, "artifact", record)
        if (isTRUE(register_run_artifacts)) {
            projectRegistryWrite(session, "run_artifact", list(
                workflow_id = identity$workflow_id,
                run_id = run_id,
                artifact_id = refs[[index]]$artifact_id,
                direction = "output",
                artifact_role = refs[[index]]$logical_key$state_role,
                ordinal = index - 1L,
                recorded_at = timestamp
            ))
        }
    }
    invisible(TRUE)
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
    identity <- context$getIdentity()
    store <- newArtifactStore(context$getPaths(), identity$project_id)
    artifactWorkflowStateEnsureMetadata(
        store,
        identity,
        protDiaArtifactDescriptorContract()
    )
    registry <- projectRegistryForContext(context, resource_policy)
    session <- initializeProjectRegistry(registry)
    tryCatch(
        artifactWorkflowStateEnsureWorkflow(session, identity),
        error = \(error) {
            closeProjectRegistry(session)
            stop(error)
        }
    )
    list(
        stage_id = stage_id,
        identity = identity,
        store = store,
        session = session,
        run_id = artifactOpaqueId("run"),
        action_id = artifactOpaqueId("action"),
        generation_id = artifactOpaqueId("gen")
    )
}

protDiaArtifactWriteStageRefs <- function(stage, tables, failure_injector = NULL) {
    tables <- lapply(names(tables), \(role) {
        protDiaArtifactPortableTable(
            tables[[role]],
            paste(stage$stage_id, role, sep = ".")
        )
    }) |>
        stats::setNames(names(tables))
    refs <- lapply(names(tables), \(role) {
        encoded <- encodeArtifactTable(
            tables[[role]],
            owner = paste("proteomics", "diann", stage$stage_id, role, sep = ".")
        )
        artifactStoreWriteParquet(
            stage$store,
            encoded,
            logical_key = list(
                project_id = stage$identity$project_id,
                omic_type = stage$identity$omic_type,
                workflow_slug = stage$identity$workflow_slug,
                stage_id = stage$stage_id,
                state_role = role,
                generation_id = stage$generation_id
            ),
            provenance_ids = c(stage$run_id, stage$action_id),
            failure_injector = failure_injector
        )
    })
    names(refs) <- names(tables)
    for (role in names(refs)) {
        protDiaArtifactVerifyRef(stage$store, refs[[role]], tables[[role]])
    }
    artifactStoreInvokeFailure(
        failure_injector,
        "before_registry_commit",
        list(stage_id = stage$stage_id, run_id = stage$run_id, refs = refs)
    )
    refs
}

protDiaArtifactRegisterStage <- function(
    stage,
    refs,
    parameters,
    source,
    deferred_commit
) {
    timestamp <- artifactRefUtcNow()
    run_status <- if (isTRUE(deferred_commit)) "running" else "completed"
    artifact_status <- if (isTRUE(deferred_commit)) "validated" else "committed"
    artifactWorkflowStateTransaction(stage$session, \() {
        protDiaArtifactRegistryRun(
            stage$session,
            stage$identity,
            stage$run_id,
            stage$action_id,
            run_status,
            timestamp
        )
        protDiaArtifactRegistrySource(
            stage$session,
            stage$identity,
            stage$run_id,
            source,
            timestamp
        )
        protDiaArtifactRegistryParameters(
            stage$session,
            stage$identity,
            stage$run_id,
            parameters,
            timestamp
        )
        protDiaArtifactRegistrySoftware(
            stage$session,
            stage$identity,
            stage$run_id,
            timestamp
        )
        protDiaArtifactRegistryRefs(
            stage$session,
            stage$identity,
            stage$run_id,
            refs,
            artifact_status,
            timestamp,
            register_run_artifacts = !isTRUE(deferred_commit)
        )
    })
    invisible(TRUE)
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
    stage <- protDiaArtifactStageResources(context, stage_id, resource_policy)
    on.exit(closeProjectRegistry(stage$session), add = TRUE)
    refs <- protDiaArtifactWriteStageRefs(stage, tables, failure_injector)
    protDiaArtifactRegisterStage(
        stage,
        refs,
        parameters,
        source,
        deferred_commit
    )
    list(
        stage_id = stage_id,
        run_id = stage$run_id,
        action_id = stage$action_id,
        generation_id = stage$generation_id,
        refs = lapply(refs, unclass),
        committed = !isTRUE(deferred_commit)
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
    log_warn = logger::log_warn
) {
    runProtDiaArtifactSafely(
        workflow_data,
        "import",
        \() {
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
