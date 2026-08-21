# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

protDiaEnrichTableDigest <- function(table) {
    artifactSemanticDigest(table)
}

protDiaEnrichFlattenTable <- function(table) {
    if (!is.data.frame(table)) {
        protDiaEnrichArtifactAbort(
            "enrichment artifact payload must be a data.frame",
            "multischolar_invalid_prot_dia_enrichment_table"
        )
    }
    if (.PROT_DIA_ENRICH_ROW_ORDER_COLUMN %in% names(table)) {
        protDiaEnrichArtifactAbort(
            "enrichment table collides with the reserved row-order column",
            "multischolar_invalid_prot_dia_enrichment_table"
        )
    }
    list_columns <- names(table)[vapply(table, is.list, logical(1))]
    list_column_classes <- lapply(list_columns, function(column) {
        unname(class(table[[column]]))
    })
    flattened <- table
    for (column in list_columns) {
        flattened[[column]] <- vapply(table[[column]], function(value) {
            jsonlite::toJSON(
                value,
                auto_unbox = TRUE,
                null = "null",
                na = "string",
                digits = NA
            )
        }, character(1))
    }
    flattened[[.PROT_DIA_ENRICH_ROW_ORDER_COLUMN]] <- seq_len(nrow(flattened))
    list(
        table = flattened,
        list_columns = as.list(list_columns),
        list_column_classes = list_column_classes,
        semantic_digest = protDiaEnrichTableDigest(list(
            table = flattened,
            list_columns = as.list(list_columns),
            list_column_classes = list_column_classes
        ))
    )
}

protDiaEnrichNormalizeListColumnClasses <- function(value, list_columns) {
    if (!is.list(value) || length(value) != length(list_columns)) {
        protDiaEnrichArtifactAbort(
            "enrichment list-column metadata is malformed",
            "multischolar_invalid_prot_dia_enrichment_manifest"
        )
    }
    lapply(value, function(classes) {
        if (is.list(classes)) classes <- unlist(classes, use.names = FALSE)
        if (!is.character(classes) || length(classes) == 0L || anyNA(classes)) {
            protDiaEnrichArtifactAbort(
                "enrichment list-column classes are malformed",
                "multischolar_invalid_prot_dia_enrichment_manifest"
            )
        }
        classes
    })
}

protDiaEnrichRestoreTable <- function(
    table,
    list_columns,
    list_column_classes
) {
    list_columns <- unlist(list_columns, use.names = FALSE)
    list_column_classes <- protDiaEnrichNormalizeListColumnClasses(
        list_column_classes,
        list_columns
    )
    for (index in seq_along(list_columns)) {
        column <- list_columns[[index]]
        if (!column %in% names(table)) next
        restored <- lapply(table[[column]], function(value) {
            restored <- jsonlite::fromJSON(value, simplifyVector = TRUE)
            if (is.list(restored) && length(restored) == 0L) character() else restored
        })
        classes <- list_column_classes[[index]]
        if (!identical(classes, "list")) class(restored) <- classes
        table[[column]] <- restored
    }
    table[[.PROT_DIA_ENRICH_ROW_ORDER_COLUMN]] <- NULL
    table
}

protDiaEnrichTableRole <- function(role, token) {
    role <- gsub("[^A-Za-z0-9_]", "_", role)
    token <- gsub("[^A-Za-z0-9_]", "_", token)
    paste0("enrichment_", role, "_", token)
}

protDiaEnrichWriteTable <- function(
    store,
    table,
    role,
    generation_id,
    token,
    provenance_ids = generation_id,
    stable_key = NULL,
    failure_injector = NULL
) {
    flattened <- protDiaEnrichFlattenTable(table)
    encoded <- encodeArtifactTable(
        flattened$table,
        stable_key = stable_key,
        owner = paste("DIA-NN enrichment", role, token)
    )
    labels <- store$labels
    logical_key <- list(
        project_id = store$project_id,
        omic_type = labels$omic_type,
        workflow_slug = labels$workflow_slug,
        stage_id = "enrichment",
        state_role = protDiaEnrichTableRole(role, token),
        generation_id = generation_id
    )
    ref <- tryCatch(
        artifactStoreWriteParquet(
            store,
            encoded,
            logical_key,
            provenance_ids = provenance_ids,
            failure_injector = failure_injector
        ),
        multischolar_duplicate_artifact_generation = \(error) {
            existing <- protDiaEnrichExistingTableRef(
                store,
                logical_key,
                encoded,
                provenance_ids
            )
            if (is.null(existing)) stop(error)
            existing
        }
    )
    list(
        role = role,
        ref = unclass(ref),
        parquet_semantic_digest = encoded$metadata$semantic_digest,
        semantic_digest = flattened$semantic_digest,
        list_columns = flattened$list_columns,
        list_column_classes = flattened$list_column_classes,
        rows = as.integer(nrow(table)),
        columns = as.integer(ncol(flattened$table))
    )
}

protDiaEnrichValidateTable <- function(entry, store) {
    required <- c(
        "role", "ref", "parquet_semantic_digest", "semantic_digest",
        "list_columns", "list_column_classes", "rows", "columns"
    )
    valid <- is.list(entry) && identical(names(entry), required) &&
        workflowCapabilityScalarString(entry$role) &&
        is.list(entry$list_columns) && is.list(entry$list_column_classes)
    if (!isTRUE(valid)) {
        protDiaEnrichArtifactAbort(
            "enrichment table reference is malformed",
            "multischolar_invalid_prot_dia_enrichment_manifest"
        )
    }
    ref <- artifactStoreNormalizeRef(entry$ref)
    sidecar <- artifactStoreReadSidecar(
        store,
        artifactStoreManagedPaths(
            store,
            ref$logical_key,
            ref$artifact_id
        )$sidecar,
        validate_payload = TRUE
    )
    valid <- identical(artifactStoreNormalizeRef(sidecar$artifact_ref), ref) &&
        identical(
            ref$hash_policy$semantic$digest,
            entry$parquet_semantic_digest
        ) && identical(as.integer(entry$rows), ref$shape$rows) &&
        identical(as.integer(entry$columns), ref$shape$columns)
    if (!isTRUE(valid)) {
        protDiaEnrichArtifactAbort(
            "enrichment table reference differs from its immutable payload",
            "multischolar_prot_dia_enrichment_table_mismatch"
        )
    }
    list_columns <- unlist(entry$list_columns, use.names = FALSE)
    entry$list_column_classes <- protDiaEnrichNormalizeListColumnClasses(
        entry$list_column_classes,
        list_columns
    )
    table <- protDiaArtifactReadTable(store, ref)
    observed_digest <- protDiaEnrichTableDigest(list(
        table = table,
        list_columns = entry$list_columns,
        list_column_classes = entry$list_column_classes
    ))
    if (!identical(observed_digest, entry$semantic_digest)) {
        protDiaEnrichArtifactAbort(
            "enrichment table semantic fingerprint differs",
            "multischolar_prot_dia_enrichment_table_mismatch"
        )
    }
    entry$ref <- unclass(ref)
    entry$rows <- as.integer(entry$rows)
    entry$columns <- as.integer(entry$columns)
    entry
}

protDiaEnrichReadTable <- function(store, entry) {
    entry <- protDiaEnrichValidateTable(entry, store)
    table <- protDiaArtifactReadTable(store, entry$ref)
    protDiaEnrichRestoreTable(
        table,
        entry$list_columns,
        entry$list_column_classes
    )
}

protDiaEnrichRequestTable <- function(request) {
    identifiers <- data.frame(
        role = "query",
        sequence = seq_along(request$identifiers),
        identifier = request$identifiers,
        stringsAsFactors = FALSE
    )
    background <- data.frame(
        role = "background",
        sequence = seq_along(request$background),
        identifier = request$background,
        stringsAsFactors = FALSE
    )
    rbind(identifiers, background)
}

protDiaEnrichRequestMetadata <- function(request) {
    list(
        backend = request$backend,
        contrast = request$contrast,
        direction = request$direction,
        query_count = as.integer(length(request$identifiers)),
        query_digest = artifactSemanticDigest(request$identifiers),
        background_count = as.integer(length(request$background)),
        background_digest = artifactSemanticDigest(request$background),
        parameters = request$parameters,
        mapping = request$mapping
    )
}

protDiaEnrichValidateResponseManifest <- function(manifest, context) {
    required <- c(
        "schema", "schema_version", "project_id", "workflow_id",
        "request_id", "request_digest", "backend", "request",
        "identifiers", "response", "response_meta", "software",
        "service_version", "created_at", "manifest_digest"
    )
    valid <- is.list(manifest) && identical(names(manifest), required) &&
        identical(manifest$schema, .PROT_DIA_ENRICH_RESPONSE_SCHEMA) &&
        identical(
            workflowStateVersionValue(manifest$schema_version),
            .PROT_DIA_ENRICH_RESPONSE_VERSION
        ) && all(vapply(
            manifest[c("project_id", "workflow_id", "request_id", "backend")],
            workflowCapabilityScalarString,
            logical(1)
        )) && is.list(manifest$request) && is.list(manifest$response_meta) &&
        is.list(manifest$software) && artifactRefValidUtc(manifest$created_at)
    if (!isTRUE(valid)) {
        protDiaEnrichArtifactAbort(
            "enrichment service response manifest is malformed",
            "multischolar_invalid_prot_dia_enrichment_manifest"
        )
    }
    if (!manifest$backend %in% c("gprofiler2", "clusterprofiler") ||
        !workflowCapabilityScalarString(manifest$service_version)) {
        protDiaEnrichArtifactAbort(
            "enrichment service response backend is unsupported",
            "multischolar_invalid_prot_dia_enrichment_manifest"
        )
    }
    artifactValidatePathComponent(manifest$request_id, "enrichment_request_id")
    artifactRefValidateDigest(
        manifest$request_digest,
        "enrichment request digest"
    )
    protDiaEnrichAssertSafeMetadata(
        manifest,
        "DIA-NN enrichment service response manifest"
    )
    if (!identical(manifest$manifest_digest, protDiaEnrichJsonDigest(manifest))) {
        protDiaEnrichArtifactAbort(
            "enrichment service response manifest fingerprint differs",
            "multischolar_prot_dia_enrichment_manifest_digest_mismatch"
        )
    }
    identity <- context$getIdentity()
    if (!identical(manifest$project_id, identity$project_id) ||
        !identical(manifest$workflow_id, identity$workflow_id)) {
        protDiaEnrichArtifactAbort(
            "enrichment response belongs to another workflow",
            "multischolar_prot_dia_enrichment_identity_mismatch"
        )
    }
    store <- protDiaEnrichStore(context)
    manifest$identifiers <- protDiaEnrichValidateTable(
        manifest$identifiers,
        store
    )
    manifest$response <- protDiaEnrichValidateTable(manifest$response, store)
    identifiers <- protDiaEnrichReadTable(store, manifest$identifiers)
    manifest$request <- protDiaEnrichValidateRequestMetadata(
        manifest$request,
        manifest$backend
    )
    request <- protDiaEnrichRequestFromManifest(manifest$request, identifiers)
    if (!identical(manifest$request_id, protDiaEnrichRequestId(request)) ||
        !identical(
            manifest$request_digest,
            artifactSemanticDigest(request)
        )) {
        protDiaEnrichArtifactAbort(
            "enrichment response request identity differs",
            "multischolar_prot_dia_enrichment_request_mismatch"
        )
    }
    manifest$response_meta <- protDiaEnrichValidateResponseMeta(
        manifest$response_meta,
        request
    )
    manifest$software <- protDiaEnrichValidateSoftware(
        manifest$software,
        manifest$backend
    )
    manifest$schema_version <- .PROT_DIA_ENRICH_RESPONSE_VERSION
    manifest
}

protDiaEnrichServiceVersion <- function(response, backend) {
    if (identical(backend, "clusterprofiler")) {
        return(protDiaNormPackageVersion("clusterProfiler"))
    }
    version <- tryCatch(response$meta$version, error = function(error) NULL)
    if (!workflowCapabilityScalarString(version)) "unreported" else version
}

protDiaEnrichPersistServiceResponse <- function(
    workflow_data,
    request,
    response,
    failure_injector = NULL,
    now = Sys.time()
) {
    context <- workflow_data$workflow_context
    request_id <- protDiaEnrichRequestId(request)
    paths <- protDiaEnrichPaths(context, request_id = request_id)
    manifest_path <- artifactResolveContainedPath(
        context$getProjectRoot(),
        paths$response_manifest
    )
    if (file.exists(manifest_path)) {
        manifest <- protDiaEnrichReadJson(
            manifest_path,
            function(value) protDiaEnrichValidateResponseManifest(value, context)
        )
        if (!identical(manifest$request_digest, artifactSemanticDigest(request))) {
            protDiaEnrichArtifactAbort(
                "cached enrichment request fingerprint differs",
                "multischolar_prot_dia_enrichment_request_mismatch"
            )
        }
        response_table <- protDiaEnrichResponseTable(request$backend, response)
        candidate <- protDiaEnrichFlattenTable(response_table)
        candidate_meta <- protDiaEnrichSafeResponseMeta(response, request)
        candidate_software <- protDiaEnrichSoftware(request$backend)
        candidate_service_version <- protDiaEnrichServiceVersion(
            response,
            request$backend
        )
        same_response <- identical(
            manifest$response$semantic_digest,
            candidate$semantic_digest
        ) && identical(
            artifactSemanticDigest(manifest$response_meta),
            artifactSemanticDigest(candidate_meta)
        ) && identical(
            artifactSemanticDigest(manifest$software),
            artifactSemanticDigest(candidate_software)
        ) && identical(
            manifest$service_version,
            candidate_service_version
        )
        if (!isTRUE(same_response)) {
            protDiaEnrichArtifactAbort(
                "immutable enrichment response differs for the same request",
                "multischolar_prot_dia_enrichment_immutable_conflict"
            )
        }
        return(list(
            request_id = request_id,
            response_digest = manifest$response$semantic_digest,
            rows = manifest$response$rows,
            manifest_relative_path = paths$response_manifest,
            manifest_digest = manifest$manifest_digest,
            result = response
        ))
    }
    store <- protDiaEnrichStore(context)
    identifiers <- protDiaEnrichWriteTable(
        store,
        protDiaEnrichRequestTable(request),
        role = "request_ids",
        generation_id = request_id,
        token = request_id,
        stable_key = c("role", "sequence"),
        failure_injector = failure_injector
    )
    response_table <- protDiaEnrichResponseTable(request$backend, response)
    result <- protDiaEnrichWriteTable(
        store,
        response_table,
        role = "service_response",
        generation_id = request_id,
        token = request_id,
        failure_injector = failure_injector
    )
    identity <- context$getIdentity()
    response_meta <- protDiaEnrichSafeResponseMeta(response, request)
    manifest <- list(
        schema = .PROT_DIA_ENRICH_RESPONSE_SCHEMA,
        schema_version = .PROT_DIA_ENRICH_RESPONSE_VERSION,
        project_id = identity$project_id,
        workflow_id = identity$workflow_id,
        request_id = request_id,
        request_digest = artifactSemanticDigest(request),
        backend = request$backend,
        request = protDiaEnrichRequestMetadata(request),
        identifiers = identifiers,
        response = result,
        response_meta = response_meta,
        software = protDiaEnrichSoftware(request$backend),
        service_version = protDiaEnrichServiceVersion(response, request$backend),
        created_at = format(now, "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC"),
        manifest_digest = NULL
    )
    manifest$manifest_digest <- protDiaEnrichJsonDigest(manifest)
    protDiaEnrichWriteJson(
        manifest,
        manifest_path,
        function(value) protDiaEnrichValidateResponseManifest(value, context),
        failure_injector = failure_injector,
        failure_stage = "before_enrichment_response_publication"
    )
    manifest <- protDiaEnrichValidateResponseManifest(manifest, context)
    list(
        request_id = request_id,
        response_digest = manifest$response$semantic_digest,
        rows = manifest$response$rows,
        manifest_relative_path = paths$response_manifest,
        manifest_digest = manifest$manifest_digest,
        result = response
    )
}

protDiaEnrichLoadServiceResponse <- function(
    workflow_data,
    request,
    required = FALSE
) {
    context <- workflow_data$workflow_context
    request_id <- protDiaEnrichRequestId(request)
    paths <- protDiaEnrichPaths(context, request_id = request_id)
    path <- artifactResolveContainedPath(
        context$getProjectRoot(), paths$response_manifest
    )
    if (!file.exists(path)) {
        if (isTRUE(required)) {
            protDiaEnrichArtifactAbort(
                "required enrichment response manifest is absent",
                "multischolar_missing_prot_dia_enrichment_replay"
            )
        }
        return(NULL)
    }
    manifest <- protDiaEnrichReadJson(
        path,
        function(value) protDiaEnrichValidateResponseManifest(value, context)
    )
    if (!identical(manifest$request_digest, artifactSemanticDigest(request))) {
        protDiaEnrichArtifactAbort(
            "replayed enrichment request does not match its manifest",
            "multischolar_prot_dia_enrichment_request_mismatch"
        )
    }
    store <- protDiaEnrichStore(context)
    ids <- protDiaEnrichReadTable(store, manifest$identifiers)
    observed_query <- ids$identifier[ids$role == "query"]
    observed_background <- ids$identifier[ids$role == "background"]
    if (!identical(observed_query, request$identifiers) ||
        !identical(observed_background, request$background)) {
        protDiaEnrichArtifactAbort(
            "replayed enrichment identifiers differ from the exact request",
            "multischolar_prot_dia_enrichment_request_mismatch"
        )
    }
    table <- protDiaEnrichReadTable(store, manifest$response)
    result <- protDiaEnrichBuildResponse(manifest, table)
    list(
        result = result,
        record = list(
            request_id = request_id,
            request_digest = manifest$request_digest,
            backend = request$backend,
            contrast = request$contrast,
            direction = request$direction,
            status = "succeeded",
            execution_state = "replay",
            attempts = 0L,
            response = list(
                request_id = request_id,
                response_digest = manifest$response$semantic_digest,
                rows = manifest$response$rows,
                manifest_relative_path = paths$response_manifest,
                manifest_digest = manifest$manifest_digest
            )
        )
    )
}

protDiaEnrichRunSource <- function(index, entry) {
    list(
        da_run_id = index$run_id,
        da_manifest_relative_path = index$manifest_relative_path,
        da_manifest_digest = index$manifest_digest,
        source_generation_id = index$source_generation_id,
        contrast_id = entry$contrast_id,
        contrast_name = entry$contrast_name,
        contrast_manifest_relative_path = entry$manifest_relative_path,
        contrast_manifest_digest = entry$manifest_digest
    )
}

protDiaEnrichValidateRecord <- function(
    record,
    context = NULL,
    parameters = NULL,
    source = NULL
) {
    required <- c(
        "request_id", "request_digest", "backend", "contrast", "direction",
        "status", "execution_state", "attempts", "response"
    )
    valid <- is.list(record) && identical(names(record), required) &&
        all(vapply(
            record[c("request_id", "backend", "contrast", "direction", "status")],
            workflowCapabilityScalarString,
            logical(1)
        )) && record$backend %in% c("gprofiler2", "clusterprofiler") &&
        record$direction %in% c("up", "down") &&
        record$status %in% c("succeeded", "skipped_empty_input") &&
        length(record$attempts) == 1L && is.numeric(record$attempts) &&
        !is.na(record$attempts) && is.finite(record$attempts) &&
        record$attempts >= 0L && record$attempts == floor(record$attempts)
    if (!isTRUE(valid)) {
        protDiaEnrichArtifactAbort(
            "enrichment execution record is incomplete or failed",
            "multischolar_invalid_prot_dia_enrichment_manifest"
        )
    }
    artifactValidatePathComponent(record$request_id, "enrichment_request_id")
    artifactRefValidateDigest(record$request_digest, "enrichment request digest")
    record$attempts <- as.integer(record$attempts)
    if (identical(record$status, "skipped_empty_input")) {
        valid <- identical(record$execution_state, "not_called") &&
            identical(record$attempts, 0L) && is.null(record$response)
        if (!isTRUE(valid)) {
            protDiaEnrichArtifactAbort(
                "skipped enrichment request has execution or response data",
                "multischolar_invalid_prot_dia_enrichment_manifest"
            )
        }
        return(record)
    }
    allowed_states <- if (identical(record$backend, "gprofiler2")) {
        c("live", "cache", "replay")
    } else {
        "local"
    }
    response_required <- c(
        "request_id", "response_digest", "rows", "manifest_relative_path",
        "manifest_digest"
    )
    response <- record$response
    valid <- record$execution_state %in% allowed_states &&
        is.list(response) && identical(names(response), response_required) &&
        identical(response$request_id, record$request_id) &&
        length(response$rows) == 1L && is.numeric(response$rows) &&
        !is.na(response$rows) && is.finite(response$rows) &&
        response$rows >= 0L && response$rows == floor(response$rows)
    if (!isTRUE(valid)) {
        protDiaEnrichArtifactAbort(
            "successful enrichment request has incomplete response provenance",
            "multischolar_invalid_prot_dia_enrichment_manifest"
        )
    }
    artifactRefValidateDigest(
        response$response_digest,
        "enrichment response digest"
    )
    artifactRefValidateDigest(
        response$manifest_digest,
        "enrichment response manifest digest"
    )
    response$rows <- as.integer(response$rows)
    response$manifest_relative_path <- artifactNormalizeRelativePath(
        response$manifest_relative_path
    )
    record$response <- response
    if (!is.null(context)) {
        expected_path <- protDiaEnrichPaths(
            context,
            request_id = record$request_id
        )$response_manifest
        if (!identical(response$manifest_relative_path, expected_path)) {
            protDiaEnrichArtifactAbort(
                "enrichment response manifest path differs from its request",
                "multischolar_prot_dia_enrichment_response_mismatch"
            )
        }
        path <- artifactResolveContainedPath(
            context$getProjectRoot(),
            response$manifest_relative_path,
            must_exist = TRUE
        )
        manifest <- protDiaEnrichReadJson(
            path,
            function(value) {
                protDiaEnrichValidateResponseManifest(value, context)
            }
        )
        valid <- identical(manifest$request_id, record$request_id) &&
            identical(manifest$request_digest, record$request_digest) &&
            identical(manifest$backend, record$backend) &&
            identical(manifest$request$contrast, record$contrast) &&
            identical(manifest$request$direction, record$direction) &&
            identical(
                manifest$response$semantic_digest,
                response$response_digest
            ) && identical(manifest$response$rows, response$rows) &&
            identical(manifest$manifest_digest, response$manifest_digest)
        if (!is.null(parameters) && !is.null(source)) {
            request_parameters <- manifest$request$parameters
            valid <- valid && identical(record$backend, parameters$backend) &&
                identical(record$contrast, source$contrast_name) &&
                identical(
                    request_parameters$user_threshold,
                    parameters$q_cutoff
                )
            if (identical(record$backend, "gprofiler2")) {
                organism <- resolveEnrichmentOrganism(
                    parameters$organism_taxid
                )
                valid <- valid && isTRUE(organism$supported) &&
                    identical(request_parameters$organism, organism$species) &&
                    identical(
                        request_parameters$correction_method,
                        parameters$correction_method
                    )
            } else {
                valid <- valid && identical(
                    request_parameters$organism_taxid,
                    parameters$organism_taxid
                )
            }
        }
        if (!isTRUE(valid)) {
            protDiaEnrichArtifactAbort(
                "enrichment response differs from its execution record",
                "multischolar_prot_dia_enrichment_response_mismatch"
            )
        }
    }
    record
}

protDiaEnrichValidateRunManifest <- function(manifest, context) {
    required <- c(
        "schema", "schema_version", "project_id", "workflow_id", "run_id",
        "source", "parameters", "parameters_digest", "software", "requests",
        "result_table", "products", "created_at", "manifest_digest"
    )
    valid <- is.list(manifest) && identical(names(manifest), required) &&
        identical(manifest$schema, .PROT_DIA_ENRICH_RUN_SCHEMA) &&
        identical(
            workflowStateVersionValue(manifest$schema_version),
            .PROT_DIA_ENRICH_RUN_VERSION
        ) && all(vapply(
            manifest[c("project_id", "workflow_id", "run_id")],
            workflowCapabilityScalarString,
            logical(1)
        )) && is.list(manifest$source) && is.list(manifest$parameters) &&
        is.list(manifest$software) && is.list(manifest$requests) &&
        is.list(manifest$products) && artifactRefValidUtc(manifest$created_at)
    if (!isTRUE(valid)) {
        protDiaEnrichArtifactAbort(
            "enrichment run manifest is malformed",
            "multischolar_invalid_prot_dia_enrichment_manifest"
        )
    }
    artifactValidatePathComponent(manifest$run_id, "enrichment_run_id")
    artifactRefValidateDigest(
        manifest$parameters_digest,
        "enrichment parameters digest"
    )
    protDiaEnrichAssertSafeMetadata(
        manifest,
        "DIA-NN enrichment run manifest"
    )
    if (!identical(manifest$manifest_digest, protDiaEnrichJsonDigest(manifest)) ||
        !identical(
            manifest$parameters_digest,
            artifactSemanticDigest(manifest$parameters)
        )) {
        protDiaEnrichArtifactAbort(
            "enrichment run manifest fingerprint differs",
            "multischolar_prot_dia_enrichment_manifest_digest_mismatch"
        )
    }
    identity <- context$getIdentity()
    if (!identical(manifest$project_id, identity$project_id) ||
        !identical(manifest$workflow_id, identity$workflow_id)) {
        protDiaEnrichArtifactAbort(
            "enrichment run belongs to another workflow",
            "multischolar_prot_dia_enrichment_identity_mismatch"
        )
    }
    source <- protDiaEnrichValidateRunSource(manifest$source, context)
    parameters <- protDiaEnrichValidateParameters(manifest$parameters)
    software <- protDiaEnrichValidateSoftware(
        manifest$software,
        parameters$backend
    )
    records <- lapply(manifest$requests, function(record) {
        protDiaEnrichValidateRecord(
            record,
            context = context,
            parameters = parameters,
            source = source
        )
    })
    records <- records[order(match(
        vapply(records, `[[`, character(1), "direction"),
        c("up", "down")
    ))]
    valid <- length(records) == 2L && identical(
        vapply(records, `[[`, character(1), "direction"),
        c("up", "down")
    ) && all(vapply(records, function(record) {
        identical(record$backend, parameters$backend) &&
            identical(record$contrast, source$contrast_name)
    }, logical(1)))
    if (!isTRUE(valid)) {
        protDiaEnrichArtifactAbort(
            "enrichment run request set is incomplete or mismatched",
            "multischolar_invalid_prot_dia_enrichment_manifest"
        )
    }
    expected <- protDiaEnrichRunId(
        source,
        parameters,
        records,
        software,
        manifest$result_table$semantic_digest
    )
    if (!identical(manifest$run_id, expected)) {
        protDiaEnrichArtifactAbort(
            "enrichment run ID differs from its scientific inputs",
            "multischolar_prot_dia_enrichment_manifest_digest_mismatch"
        )
    }
    store <- protDiaEnrichStore(context)
    manifest$result_table <- protDiaEnrichValidateTable(
        manifest$result_table,
        store
    )
    result_table <- protDiaEnrichReadTable(store, manifest$result_table)
    protDiaEnrichValidateResultSemantics(result_table, parameters, records)
    manifest$products <- protDiaEnrichValidateProducts(
        manifest$products,
        source,
        records,
        context,
        manifest$run_id
    )
    manifest$source <- source
    manifest$parameters <- parameters
    manifest$software <- software
    manifest$requests <- records
    manifest$schema_version <- .PROT_DIA_ENRICH_RUN_VERSION
    manifest
}

protDiaEnrichPersistRun <- function(
    workflow_data,
    selected,
    parameters,
    records,
    result_table,
    pathway_dir,
    failure_injector = NULL,
    now = Sys.time()
) {
    context <- workflow_data$workflow_context
    source <- protDiaEnrichValidateRunSource(
        protDiaEnrichRunSource(selected$index, selected$entry),
        context
    )
    parameters <- protDiaEnrichValidateParameters(parameters)
    software <- protDiaEnrichValidateSoftware(
        protDiaEnrichSoftware(parameters$backend),
        parameters$backend
    )
    records <- lapply(records, function(record) {
        protDiaEnrichValidateRecord(
            record,
            context = context,
            parameters = parameters,
            source = source
        )
    })
    records <- records[order(match(
        vapply(records, `[[`, character(1), "direction"),
        c("up", "down")
    ))]
    result_digest <- protDiaEnrichFlattenTable(result_table)$semantic_digest
    run_id <- protDiaEnrichRunId(
        source,
        parameters,
        records,
        software,
        result_digest
    )
    paths <- protDiaEnrichPaths(context, run_id = run_id)
    path <- artifactResolveContainedPath(
        context$getProjectRoot(), paths$run_manifest
    )
    if (file.exists(path)) {
        manifest <- protDiaEnrichReadJson(
            path,
            function(value) protDiaEnrichValidateRunManifest(value, context)
        )
        return(list(manifest = manifest, resumed = TRUE))
    }
    store <- protDiaEnrichStore(context)
    table_entry <- protDiaEnrichWriteTable(
        store,
        result_table,
        role = "display_results",
        generation_id = run_id,
        token = selected$entry$contrast_id,
        provenance_ids = c(run_id, source$da_run_id),
        failure_injector = failure_injector
    )
    products <- protDiaEnrichPersistProducts(
        context,
        run_id,
        pathway_dir,
        selected$entry$contrast_name,
        failure_injector
    )
    identity <- context$getIdentity()
    manifest <- list(
        schema = .PROT_DIA_ENRICH_RUN_SCHEMA,
        schema_version = .PROT_DIA_ENRICH_RUN_VERSION,
        project_id = identity$project_id,
        workflow_id = identity$workflow_id,
        run_id = run_id,
        source = source,
        parameters = parameters,
        parameters_digest = artifactSemanticDigest(parameters),
        software = software,
        requests = records,
        result_table = table_entry,
        products = products,
        created_at = format(now, "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC"),
        manifest_digest = NULL
    )
    manifest$manifest_digest <- protDiaEnrichJsonDigest(manifest)
    protDiaEnrichWriteJson(
        manifest,
        path,
        function(value) protDiaEnrichValidateRunManifest(value, context),
        failure_injector = failure_injector,
        failure_stage = "before_enrichment_run_publication"
    )
    list(
        manifest = protDiaEnrichValidateRunManifest(manifest, context),
        resumed = FALSE
    )
}

publishProtDiaEnrichRun <- function(
    workflow_data,
    prepared,
    failure_injector = NULL
) {
    context <- workflow_data$workflow_context
    generation <- workflow_data$state_manager$getCurrentGenerationId()
    if (!identical(
        prepared$manifest$source$source_generation_id,
        generation
    )) {
        protDiaEnrichArtifactAbort(
            "enrichment source generation changed before publication",
            "multischolar_stale_prot_dia_enrichment_da_index"
        )
    }
    current <- artifactResolveContainedPath(
        context$getProjectRoot(), protDiaEnrichPaths(context)$current
    )
    protDiaEnrichWriteJson(
        prepared$manifest,
        current,
        function(value) protDiaEnrichValidateRunManifest(value, context),
        replace = TRUE,
        failure_injector = failure_injector,
        failure_stage = "before_enrichment_current_publication"
    )
    prepared$current_path <- current
    prepared$published <- TRUE
    prepared
}
