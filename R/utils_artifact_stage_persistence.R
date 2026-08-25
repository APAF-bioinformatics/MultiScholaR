# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

#' Abort shared artifact stage persistence
#'
#' @param message Human-readable failure message.
#' @param class Specific condition class.
#' @param ... Additional fields passed to [rlang::abort()].
#'
#' @return This function does not return; it signals a typed error.
#' @noRd
artifactStagePersistenceAbort <- function(message, class, ...) {
    rlang::abort(
        message,
        class = c(class, "multischolar_artifact_stage_persistence_error"),
        ...
    )
}

#' Extract the immutable descriptor contract
#'
#' @param descriptor Workflow descriptor to validate.
#'
#' @return Descriptor identifier, version, and digest fields.
#' @noRd
artifactStageDescriptorContract <- function(descriptor) {
    descriptor <- validateArtifactWorkflowDescriptor(descriptor)
    list(
        descriptor_id = descriptor$descriptor_id,
        descriptor_version = descriptor$descriptor_version,
        descriptor_digest = descriptor$descriptor_digest
    )
}

#' Test an identity against an exact workflow descriptor
#'
#' @param identity Candidate workflow identity.
#' @param descriptor Workflow descriptor to validate.
#'
#' @return A scalar logical indicating an exact identity match.
#' @noRd
artifactStageIdentityMatches <- function(identity, descriptor) {
    if (is.null(identity)) return(FALSE)
    descriptor <- validateArtifactWorkflowDescriptor(descriptor)
    identical(
        workflowCapabilityKey(identity),
        workflowCapabilityKey(descriptor$identity)
    )
}

#' Test whether the coordinator owns an exact artifact context
#'
#' @param workflow_data Mutable workflow state.
#' @param descriptor Workflow descriptor to validate.
#'
#' @return A scalar logical indicating artifact ownership and identity parity.
#' @noRd
artifactStageCoordinatorOwned <- function(workflow_data, descriptor) {
    context <- workflow_data$workflow_context
    if (!inherits(context, "WorkflowContext") || !context$isBound()) {
        return(FALSE)
    }
    decision <- context$getStorageDecision()
    identical(decision$effective_backend, "artifact") &&
        identical(decision$capability_id, descriptor$descriptor_id) &&
        identical(decision$capability_version, descriptor$descriptor_version) &&
        artifactStageIdentityMatches(context$getIdentity(), descriptor)
}

#' Resolve one exact workflow descriptor for an import tuple
#'
#' @param workflow_type Workflow type identifier.
#' @param input_format Input format identifier.
#' @param data_level Scientific data level.
#' @param descriptor_catalogue Workflow descriptor catalogue.
#'
#' @return The matching descriptor, or `NULL` when no tuple is certified.
#' @noRd
artifactStageDescriptorForImport <- function(
    workflow_type,
    input_format,
    data_level,
    descriptor_catalogue = artifactWorkflowDescriptorCatalogue()
) {
    descriptors <- artifactDescriptorCatalogueValues(descriptor_catalogue)
    matches <- vapply(descriptors, \(descriptor) {
        identity <- descriptor$identity
        identical(identity$workflow_type, workflow_type) &&
            identical(identity$input_format, input_format) &&
            identical(identity$data_level, data_level)
    }, logical(1))
    if (sum(matches) > 1L) {
        artifactStagePersistenceAbort(
            "artifact import identity matches multiple workflow descriptors",
            "multischolar_ambiguous_artifact_stage_descriptor"
        )
    }
    if (!any(matches)) return(NULL)
    descriptors[[which(matches)]]
}

#' Prepare a descriptor-bound artifact stage context
#'
#' @param workflow_data Mutable workflow state.
#' @param workflow_type Workflow type identifier.
#' @param input_format Input format identifier.
#' @param data_level Scientific data level.
#' @param descriptor_catalogue Workflow descriptor catalogue.
#'
#' @return A prepared context result with an explicit enabled flag and reason.
#' @noRd
prepareArtifactStageContext <- function(
    workflow_data,
    workflow_type,
    input_format,
    data_level,
    descriptor_catalogue = artifactWorkflowDescriptorCatalogue()
) {
    descriptor <- artifactStageDescriptorForImport(
        workflow_type,
        input_format,
        data_level,
        descriptor_catalogue
    )
    if (is.null(descriptor)) {
        return(list(enabled = FALSE, reason = "artifact_descriptor_unavailable"))
    }
    context <- workflow_data$workflow_context
    if (!inherits(context, "WorkflowContext")) {
        return(list(enabled = FALSE, reason = "workflow_context_absent"))
    }
    if (!context$isBound()) {
        bindWorkflowContextFromImport(
            context,
            workflow_type = workflow_type,
            input_format = input_format,
            data_level = data_level,
            descriptor_catalogue = descriptor_catalogue
        )
    }
    decision <- context$getStorageDecision()
    enabled <- identical(decision$effective_backend, "artifact") &&
        decision$effective_rollout %in% .WORKFLOW_ARTIFACT_ROLLOUTS &&
        identical(decision$capability_id, descriptor$descriptor_id) &&
        identical(decision$capability_version, descriptor$descriptor_version) &&
        artifactStageIdentityMatches(context$getIdentity(), descriptor)
    list(
        enabled = enabled,
        reason = decision$reason_code,
        context = context,
        decision = unclass(decision),
        descriptor = descriptor
    )
}

#' Record an additive artifact stage result
#'
#' @param workflow_data Mutable workflow state.
#' @param stage_id Stage identifier.
#' @param result Stage operation result.
#'
#' @return `result`, invisibly.
#' @noRd
recordArtifactStageResult <- function(workflow_data, stage_id, result) {
    results <- workflow_data$artifact_stage_results
    if (is.null(results)) results <- list()
    results[[stage_id]] <- result
    workflow_data$artifact_stage_results <- results
    invisible(result)
}

#' Run an additive artifact stage without replacing memory state
#'
#' @param workflow_data Mutable workflow state.
#' @param stage_id Stage identifier.
#' @param operation Zero-argument artifact operation.
#' @param owner_label Human-readable workflow owner label.
#' @param log_warn Warning logger used for additive artifact failures.
#'
#' @return The recorded operation result, invisibly.
#' @noRd
runArtifactStageSafely <- function(
    workflow_data,
    stage_id,
    operation,
    owner_label,
    log_warn = logger::log_warn
) {
    result <- tryCatch(
        operation(),
        error = \(error) {
            log_warn(paste(
                owner_label,
                "artifact",
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
    recordArtifactStageResult(workflow_data, stage_id, result)
}

#' Validate a portable artifact table
#'
#' @param value Candidate rectangular value.
#' @param owner Human-readable table owner.
#' @param abort_fn Typed abort function.
#'
#' @return A data frame suitable for artifact encoding.
#' @noRd
artifactStagePortableTable <- function(value, owner, abort_fn) {
    if (!is.data.frame(value)) {
        abort_fn(
            sprintf("artifact '%s' must be a data frame", owner),
            "multischolar_invalid_artifact_stage_table"
        )
    }
    if (inherits(value, "spec_tbl_df")) value <- tibble::as_tibble(value)
    value
}

#' Encode named metadata as a rectangular table
#'
#' @param value Named metadata value or list.
#'
#' @return A data frame of keys, canonical JSON values, and semantic digests.
#' @noRd
artifactStageMetadataTable <- function(value) {
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
        value_json = vapply(value, artifactStageParameterJson, character(1)),
        value_digest = vapply(value, artifactSemanticDigest, character(1)),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
}

#' Normalize an optional artifact table
#'
#' @param value Optional rectangular value.
#' @param role Artifact state role.
#' @param abort_fn Typed abort function.
#'
#' @return The input table or an explicit empty availability table.
#' @noRd
artifactStageOptionalTable <- function(value, role, abort_fn) {
    if (is.null(value)) {
        return(data.frame(
            artifact_role = character(),
            availability = character(),
            stringsAsFactors = FALSE
        ))
    }
    artifactStagePortableTable(value, role, abort_fn)
}

#' Classify a contrasts payload for provenance
#'
#' @param value Contrasts payload.
#'
#' @return One of `"null"`, `"character"`, or `"data.frame"`.
#' @noRd
artifactStageContrastsKind <- function(value) {
    if (is.null(value)) return("null")
    if (is.character(value)) return("character")
    "data.frame"
}

#' Resolve the current MultiScholaR package version
#'
#' @return The installed or source package version as a scalar character value.
#' @noRd
artifactStagePackageVersion <- function() {
    version <- tryCatch(
        as.character(utils::packageVersion("MultiScholaR")),
        error = \(...) NA_character_
    )
    if (!is.na(version)) return(version)
    description <- tryCatch(read.dcf("DESCRIPTION"), error = \(...) NULL)
    if (is.null(description)) "unknown" else description[[1L, "Version"]]
}

#' Serialize one provenance parameter deterministically
#'
#' @param value Parameter value.
#'
#' @return A scalar JSON character value.
#' @noRd
artifactStageParameterJson <- function(value) {
    as.character(jsonlite::toJSON(
        value,
        auto_unbox = TRUE,
        null = "null",
        na = "string",
        digits = NA
    ))
}

#' Independently hydrate and verify an artifact reference
#'
#' @param store Validated artifact store.
#' @param ref Artifact reference to hydrate.
#' @param expected Expected rectangular value.
#' @param abort_fn Typed abort function.
#'
#' @return The independently hydrated value, invisibly.
#' @noRd
artifactStageVerifyRef <- function(store, ref, expected, abort_fn) {
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
        abort_fn(
            sprintf(
                "independent artifact hydration changed '%s'",
                ref$logical_key$state_role
            ),
            "multischolar_inexact_artifact_stage_hydration"
        )
    }
    invisible(hydrated)
}

#' Register an artifact workflow run
#'
#' @param session Open project registry session.
#' @param identity Bound workflow identity.
#' @param run_id Opaque run identifier.
#' @param action_id Opaque action identifier.
#' @param status Initial run status.
#' @param timestamp Registry timestamp.
#'
#' @return The registry write result, invisibly.
#' @noRd
artifactStageRegistryRun <- function(
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

#' Register artifact source provenance
#'
#' @param session Open project registry session.
#' @param identity Bound workflow identity.
#' @param run_id Opaque run identifier.
#' @param source Optional source provenance record.
#' @param timestamp Registry timestamp.
#'
#' @return A scalar logical indicating whether a source was registered, invisibly.
#' @noRd
artifactStageRegistrySource <- function(
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

#' Register artifact parameter provenance
#'
#' @param session Open project registry session.
#' @param identity Bound workflow identity.
#' @param run_id Opaque run identifier.
#' @param parameters Optional named parameter list.
#' @param timestamp Registry timestamp.
#'
#' @return A scalar logical indicating whether parameters were registered, invisibly.
#' @noRd
artifactStageRegistryParameters <- function(
    session,
    identity,
    run_id,
    parameters,
    timestamp
) {
    if (is.null(parameters)) return(invisible(FALSE))
    for (parameter_name in names(parameters)) {
        value <- parameters[[parameter_name]]
        projectRegistryWrite(session, "parameter", list(
            workflow_id = identity$workflow_id,
            run_id = run_id,
            parameter_id = artifactOpaqueId("parameter"),
            parameter_name = parameter_name,
            value_json = artifactStageParameterJson(value),
            value_digest = artifactSemanticDigest(value),
            recorded_at = timestamp
        ))
    }
    invisible(TRUE)
}

#' Register MultiScholaR software provenance
#'
#' @param session Open project registry session.
#' @param identity Bound workflow identity.
#' @param run_id Opaque run identifier.
#' @param timestamp Registry timestamp.
#'
#' @return The registry write result, invisibly.
#' @noRd
artifactStageRegistrySoftware <- function(session, identity, run_id, timestamp) {
    version <- artifactStagePackageVersion()
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

#' Register artifact references and run dependencies
#'
#' @param session Open project registry session.
#' @param identity Bound workflow identity.
#' @param run_id Opaque run identifier.
#' @param refs Named artifact references.
#' @param status Initial artifact status.
#' @param timestamp Registry timestamp.
#' @param register_run_artifacts Whether to expose run-artifact dependencies now.
#'
#' @return `TRUE`, invisibly.
#' @noRd
artifactStageRegistryRefs <- function(
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

#' Open artifact store and registry resources for a stage
#'
#' @param context Exact bound workflow context.
#' @param descriptor Exact workflow descriptor.
#' @param stage_id Stage identifier.
#' @param resource_policy Optional project registry resource policy.
#' @param abort_fn Typed abort function.
#'
#' @return A stage resource bundle with store, session, and opaque identifiers.
#' @noRd
artifactStageResources <- function(
    context,
    descriptor,
    stage_id,
    resource_policy = NULL,
    abort_fn = artifactStagePersistenceAbort
) {
    descriptor <- validateArtifactWorkflowDescriptor(descriptor)
    valid_context <- inherits(context, "WorkflowContext") && context$isBound() &&
        artifactStageIdentityMatches(context$getIdentity(), descriptor)
    if (!isTRUE(valid_context)) {
        abort_fn(
            "artifact stage persistence requires its exact bound context",
            "multischolar_invalid_artifact_stage_context"
        )
    }
    identity <- context$getIdentity()
    store <- newArtifactStore(context$getPaths(), identity$project_id)
    registry_identity <- artifactWorkflowStateEnsureMetadata(
        store,
        identity,
        artifactStageDescriptorContract(descriptor)
    )
    registry <- projectRegistryForContext(context, resource_policy)
    session <- initializeProjectRegistry(registry)
    tryCatch(
        artifactWorkflowStateEnsureWorkflow(session, registry_identity),
        error = \(error) {
            closeProjectRegistry(session)
            stop(error)
        }
    )
    list(
        stage_id = stage_id,
        identity = identity,
        registry_identity = registry_identity,
        descriptor = descriptor,
        store = store,
        session = session,
        run_id = artifactOpaqueId("run"),
        action_id = artifactOpaqueId("action"),
        generation_id = artifactOpaqueId("gen")
    )
}

#' Atomically write and verify all rectangular stage artifacts
#'
#' @param stage Open artifact stage resource bundle.
#' @param tables Named rectangular stage tables.
#' @param failure_injector Optional artifact failure injector used by tests.
#' @param verification Artifact store verification mode.
#' @param abort_fn Typed abort function.
#'
#' @return Named immutable artifact references.
#' @noRd
artifactStageWriteRefs <- function(
    stage,
    tables,
    failure_injector = NULL,
    verification = c("inline_exact", "deferred_exact"),
    abort_fn = artifactStagePersistenceAbort
) {
    verification <- match.arg(verification)
    tables <- lapply(names(tables), \(role) {
        artifactStagePortableTable(
            tables[[role]],
            paste(stage$stage_id, role, sep = "."),
            abort_fn
        )
    }) |>
        stats::setNames(names(tables))
    refs <- lapply(names(tables), \(role) {
        encoded <- encodeArtifactTable(
            tables[[role]],
            owner = paste(
                stage$descriptor$descriptor_id,
                stage$stage_id,
                role,
                sep = "."
            )
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
            failure_injector = failure_injector,
            verification = verification
        )
    })
    names(refs) <- names(tables)
    if (identical(verification, "inline_exact")) {
        for (role in names(refs)) {
            artifactStageVerifyRef(stage$store, refs[[role]], tables[[role]], abort_fn)
        }
    }
    artifactStoreInvokeFailure(
        failure_injector,
        "before_registry_commit",
        list(stage_id = stage$stage_id, run_id = stage$run_id, refs = refs)
    )
    refs
}

#' Register one artifact stage generation transactionally
#'
#' @param stage Open artifact stage resource bundle.
#' @param refs Named immutable artifact references.
#' @param parameters Optional named provenance parameters.
#' @param source Optional source provenance record.
#' @param deferred_commit Whether visibility waits for an S4 checkpoint.
#'
#' @return `TRUE`, invisibly.
#' @noRd
artifactStageRegister <- function(
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
        artifactStageRegistryRun(
            stage$session,
            stage$registry_identity,
            stage$run_id,
            stage$action_id,
            run_status,
            timestamp
        )
        artifactStageRegistrySource(
            stage$session,
            stage$registry_identity,
            stage$run_id,
            source,
            timestamp
        )
        artifactStageRegistryParameters(
            stage$session,
            stage$registry_identity,
            stage$run_id,
            parameters,
            timestamp
        )
        artifactStageRegistrySoftware(
            stage$session,
            stage$registry_identity,
            stage$run_id,
            timestamp
        )
        artifactStageRegistryRefs(
            stage$session,
            stage$registry_identity,
            stage$run_id,
            refs,
            artifact_status,
            timestamp,
            register_run_artifacts = !isTRUE(deferred_commit)
        )
    })
    invisible(TRUE)
}

#' Write, verify, and register one descriptor-bound artifact stage
#'
#' @param context Exact bound workflow context.
#' @param descriptor Exact workflow descriptor.
#' @param stage_id Stage identifier.
#' @param tables Named rectangular stage tables.
#' @param parameters Optional named provenance parameters.
#' @param source Optional source provenance record.
#' @param deferred_commit Whether visibility waits for an S4 checkpoint.
#' @param failure_injector Optional artifact failure injector used by tests.
#' @param resource_policy Optional project registry resource policy.
#' @param abort_fn Typed abort function.
#'
#' @return A stage record containing immutable references and transaction IDs.
#' @noRd
writeArtifactStage <- function(
    context,
    descriptor,
    stage_id,
    tables,
    parameters,
    source = NULL,
    deferred_commit = FALSE,
    failure_injector = NULL,
    resource_policy = NULL,
    abort_fn = artifactStagePersistenceAbort
) {
    stage <- artifactStageResources(
        context,
        descriptor,
        stage_id,
        resource_policy,
        abort_fn
    )
    on.exit(closeProjectRegistry(stage$session), add = TRUE)
    refs <- artifactStageWriteRefs(
        stage,
        tables,
        failure_injector,
        abort_fn = abort_fn
    )
    artifactStageRegister(stage, refs, parameters, source, deferred_commit)
    list(
        stage_id = stage_id,
        run_id = stage$run_id,
        action_id = stage$action_id,
        generation_id = stage$generation_id,
        refs = lapply(refs, unclass),
        committed = !isTRUE(deferred_commit)
    )
}

#' Advance one deferred artifact run with compare-and-set semantics
#'
#' @param connection Writable project registry connection.
#' @param identity Bound workflow identity.
#' @param stage Deferred artifact stage record.
#' @param completed Whether the stage completed successfully.
#' @param timestamp Registry timestamp.
#' @param abort_fn Typed abort function.
#'
#' @return `TRUE`, invisibly.
#' @noRd
artifactStageSetRunStatus <- function(
    connection,
    identity,
    stage,
    completed,
    timestamp,
    abort_fn = artifactStagePersistenceAbort
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
        abort_fn(
            "artifact stage run status did not advance exactly once",
            "multischolar_artifact_stage_compare_and_set_failed"
        )
    }
    invisible(TRUE)
}

#' Commit a complete deferred artifact reference set
#'
#' @param session Open project registry session.
#' @param connection Writable project registry connection.
#' @param identity Bound workflow identity.
#' @param stage Deferred artifact stage record.
#' @param timestamp Registry timestamp.
#' @param abort_fn Typed abort function.
#'
#' @return `TRUE`, invisibly.
#' @noRd
artifactStageCommitRefs <- function(
    session,
    connection,
    identity,
    stage,
    timestamp,
    abort_fn = artifactStagePersistenceAbort
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
        abort_fn(
            "artifact stage refs did not commit as one complete set",
            "multischolar_artifact_stage_compare_and_set_failed"
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

#' Transactionally complete or fail a deferred artifact stage
#'
#' @param context Exact bound workflow context.
#' @param stage Deferred artifact stage record.
#' @param completed Whether the stage completed successfully.
#' @param failure_injector Optional artifact failure injector used by tests.
#' @param resource_policy Optional project registry resource policy.
#' @param abort_fn Typed abort function.
#'
#' @return `TRUE`, invisibly.
#' @noRd
artifactStageUpdateStatus <- function(
    context,
    stage,
    completed,
    failure_injector = NULL,
    resource_policy = NULL,
    abort_fn = artifactStagePersistenceAbort
) {
    identity <- context$getIdentity()
    store <- newArtifactStore(context$getPaths(), identity$project_id)
    registry_identity <- artifactRegistryIdentity(store, identity)
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
        artifactStageSetRunStatus(
            connection,
            registry_identity,
            stage,
            completed,
            timestamp,
            abort_fn
        )
        if (isTRUE(completed)) {
            artifactStageCommitRefs(
                session,
                connection,
                registry_identity,
                stage,
                timestamp,
                abort_fn
            )
        }
    })
    invisible(TRUE)
}
