# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

.METAB_IMPORT_PENDING_STAGE_SCHEMA <-
    "multischolar.metabolomics_pending_import_stage"

#' Resolve one custom metabolomics source delimiter
#'
#' @param source_path Existing custom tabular source.
#'
#' @return A one-character delimiter.
#' @noRd
metabImportStreamingDelimiter <- function(source_path) {
    if (identical(tolower(tools::file_ext(source_path)), "csv")) "," else "\t"
}

#' Quote one DuckDB identifier
#'
#' @param connection Open DuckDB connection.
#' @param value Identifier value.
#'
#' @return A quoted SQL identifier.
#' @noRd
metabImportStreamingQuote <- function(connection, value) {
    as.character(DBI::dbQuoteIdentifier(connection, value))
}

#' Read one bounded custom metabolomics header
#'
#' @param source_path Existing custom tabular source.
#' @param maximum_bytes Maximum accepted header bytes.
#'
#' @return Exact source column names.
#' @noRd
metabImportStreamingHeader <- function(
    source_path,
    maximum_bytes = 1024^2
) {
    connection <- file(source_path, open = "rb")
    on.exit(close(connection), add = TRUE)
    header <- readLines(connection, n = 1L, warn = FALSE)
    if (length(header) != 1L || nchar(header, type = "bytes") > maximum_bytes) {
        metabArtifactAbort(
            "custom metabolomics header is missing or exceeds its bound",
            "multischolar_invalid_metabolomics_ingress_header"
        )
    }
    delimiter <- metabImportStreamingDelimiter(source_path)
    parsed <- utils::read.table(
        text = header,
        header = FALSE,
        sep = delimiter,
        quote = "\"",
        comment.char = "",
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    columns <- as.character(parsed[1L, , drop = TRUE])
    if (!length(columns) || any(!nzchar(columns)) || anyDuplicated(columns)) {
        metabArtifactAbort(
            "custom metabolomics columns must be non-empty and unique",
            "multischolar_invalid_metabolomics_ingress_schema"
        )
    }
    columns
}

#' Validate one custom metabolomics streaming mapping
#'
#' @param source_path Existing source file.
#' @param column_mapping Established workflow column mapping.
#' @param sanitize_names Whether sample names are sanitized.
#'
#' @return Source and output schema contract.
#' @noRd
metabImportStreamingPreflight <- function(
    source_path,
    column_mapping,
    sanitize_names = FALSE
) {
    columns <- metabImportStreamingHeader(source_path)
    samples <- unname(as.character(column_mapping$sample_columns))
    required <- c(
        column_mapping$metabolite_id_col,
        column_mapping$annotation_col,
        samples
    )
    required <- required[!is.na(required) & nzchar(required)]
    if (!length(samples) || anyDuplicated(samples) ||
        length(setdiff(required, columns))) {
        metabArtifactAbort(
            "custom metabolomics mapping differs from its source schema",
            "multischolar_invalid_metabolomics_ingress_mapping"
        )
    }
    output_names <- columns
    if (isTRUE(sanitize_names)) {
        cleaned <- janitor::make_clean_names(samples)
        output_names[match(samples, columns)] <- cleaned
        samples <- cleaned
    }
    if (anyDuplicated(output_names)) {
        metabArtifactAbort(
            "custom metabolomics columns collide after sample sanitization",
            "multischolar_invalid_metabolomics_ingress_mapping"
        )
    }
    output_mapping <- column_mapping
    output_mapping$sample_columns <- samples
    list(
        source_size_bytes = unname(as.numeric(file.info(source_path)$size)),
        source_digest = artifactByteDigest(source_path),
        columns = columns,
        output_names = output_names,
        source_sample_columns = unname(as.character(
            column_mapping$sample_columns
        )),
        column_mapping = output_mapping,
        complete_payload_materialized = FALSE
    )
}

#' Build one ordered custom metabolomics streaming query
#'
#' @param connection Open DuckDB connection.
#' @param source_path Existing source file.
#' @param preflight Validated source preflight.
#'
#' @return Ordered canonical SQL.
#' @noRd
metabImportStreamingSelect <- function(connection, source_path, preflight) {
    delimiter <- metabImportStreamingDelimiter(source_path)
    source_literal <- as.character(DBI::dbQuoteString(connection, source_path))
    source_table <- metabImportStreamingQuote(
        connection,
        ".multischolar_metab_source"
    )
    source_order <- metabImportStreamingQuote(connection, ".source_row")
    DBI::dbExecute(connection, paste0(
        "CREATE TEMP TABLE ", source_table, " AS SELECT ",
        "row_number() OVER ()::BIGINT AS ", source_order,
        ", * FROM read_csv(", source_literal,
        ", delim = ", as.character(DBI::dbQuoteString(connection, delimiter)),
        ", header = true, auto_detect = true, sample_size = -1, ",
        "nullstr = ['', 'NA'], strict_mode = true)"
    ))
    projections <- vapply(seq_along(preflight$columns), \(index) {
        source_name <- preflight$columns[[index]]
        output_name <- preflight$output_names[[index]]
        source <- metabImportStreamingQuote(connection, source_name)
        target <- metabImportStreamingQuote(connection, output_name)
        expression <- if (source_name %in% preflight$source_sample_columns) {
            paste0("TRY_CAST(", source, " AS DOUBLE)")
        } else {
            source
        }
        paste(expression, "AS", target)
    }, character(1))
    paste0(
        "SELECT row_number() OVER (ORDER BY ", source_order,
        ")::BIGINT AS ", metabImportStreamingQuote(
            connection,
            .artifactRowOrderColumn
        ), ", ", paste(projections, collapse = ", "), " FROM ",
        source_table, " ORDER BY ", source_order
    )
}

#' Stream one custom metabolomics assay to canonical Parquet
#'
#' @param source_path Existing source file.
#' @param output_path Destination Parquet path.
#' @param column_mapping Established column mapping.
#' @param sanitize_names Whether sample names are sanitized.
#' @param row_group_rows Parquet row-group size.
#' @param memory_limit_bytes DuckDB memory limit.
#'
#' @return Source binding, output mapping, and output path.
#' @noRd
writeMetabImportStreamingParquet <- function(
    source_path,
    output_path,
    column_mapping,
    sanitize_names = FALSE,
    row_group_rows = 65536L,
    memory_limit_bytes = 128 * 1024^2
) {
    preflight <- metabImportStreamingPreflight(
        source_path,
        column_mapping,
        sanitize_names
    )
    database <- tempfile("metab-import-streaming-", fileext = ".duckdb")
    temporary <- tempfile("metab-import-spill-")
    dir.create(temporary)
    on.exit(
        unlink(c(database, temporary), recursive = TRUE, force = TRUE),
        add = TRUE
    )
    connection <- DBI::dbConnect(duckdb::duckdb(database, shared_home = FALSE))
    on.exit(DBI::dbDisconnect(connection, shutdown = TRUE), add = TRUE)
    DBI::dbExecute(connection, "SET threads = 1")
    DBI::dbExecute(connection, "SET preserve_insertion_order = true")
    DBI::dbExecute(
        connection,
        paste0("SET memory_limit = '", as.integer(memory_limit_bytes), "B'")
    )
    DBI::dbExecute(connection, paste0(
        "SET temp_directory = ",
        as.character(DBI::dbQuoteString(connection, temporary))
    ))
    select <- metabImportStreamingSelect(connection, source_path, preflight)
    destination <- as.character(DBI::dbQuoteString(connection, output_path))
    compression <- if (arrow::codec_is_available("zstd")) "ZSTD" else "SNAPPY"
    DBI::dbExecute(connection, paste0(
        "COPY (", select, ") TO ", destination,
        " (FORMAT PARQUET, COMPRESSION ", compression,
        ", ROW_GROUP_SIZE ", as.integer(row_group_rows), ")"
    ))
    list(
        path = output_path,
        preflight = preflight,
        column_mapping = preflight$column_mapping
    )
}

#' Resolve exact custom metabolomics source bindings
#'
#' @param assay1_file First assay source.
#' @param assay1_name First assay name.
#' @param assay2_file Optional second assay source.
#' @param assay2_name Optional second assay name.
#'
#' @return Named source paths in semantic assay order.
#' @noRd
metabImportStageSources <- function(
    assay1_file,
    assay1_name,
    assay2_file = NULL,
    assay2_name = NULL
) {
    paths <- assay1_file
    names(paths) <- assay1_name
    if (!is.null(assay2_file) && nzchar(assay2_file) &&
        !is.null(assay2_name) && nzchar(assay2_name)) {
        paths <- c(paths, stats::setNames(assay2_file, assay2_name))
    }
    valid <- length(paths) > 0L && !anyDuplicated(names(paths)) &&
        all(vapply(paths, \(path) {
            workflowCapabilityScalarString(path) && file.exists(path) &&
                !dir.exists(path)
        }, logical(1)))
    if (!isTRUE(valid)) {
        metabArtifactAbort(
            "custom metabolomics staging requires ordered regular sources",
            "multischolar_invalid_metabolomics_stage_sources"
        )
    }
    paths
}

#' Prepare one payload-free custom metabolomics import stage
#'
#' @param prepared Prepared exact workflow context.
#' @param resource_policy Optional project registry policy.
#'
#' @return A stage bundle without open registry resources.
#' @noRd
metabImportPreparedStage <- function(prepared, resource_policy = NULL) {
    stage <- artifactStageResources(
        prepared$context,
        prepared$descriptor,
        "import",
        resource_policy,
        metabArtifactAbort
    )
    closeProjectRegistry(stage$session)
    stage$session <- NULL
    stage$resource_policy <- normalizeProjectRegistryPolicy(resource_policy)
    artifactResourceDataOnly(stage, "metabolomics pending import stage")
    stage
}

#' Publish one streaming assay ref into a pending generation
#'
#' @param stage Prepared import stage.
#' @param role Exact assay role.
#' @param source_path Existing custom source.
#' @param column_mapping Established mapping.
#' @param sanitize_names Whether sample names are sanitized.
#'
#' @return Immutable streaming artifact reference.
#' @noRd
metabImportStageAssayRef <- function(
    stage,
    role,
    source_path,
    column_mapping,
    sanitize_names
) {
    output_path <- tempfile("metab-import-streaming-", fileext = ".parquet")
    on.exit(unlink(output_path, force = TRUE), add = TRUE)
    writeMetabImportStreamingParquet(
        source_path,
        output_path,
        column_mapping,
        sanitize_names
    )
    encoded <- encodeArtifactStreamingParquet(
        output_path,
        owner = paste(stage$descriptor$descriptor_id, "import", role, sep = "."),
        data_frame_class = "data.frame"
    )
    artifactStorePublishStreamingParquet(
        stage$store,
        encoded,
        logical_key = list(
            project_id = stage$identity$project_id,
            omic_type = stage$identity$omic_type,
            workflow_slug = stage$identity$workflow_slug,
            stage_id = "import",
            state_role = role,
            generation_id = stage$generation_id
        ),
        provenance_ids = c(stage$run_id, stage$action_id),
        verification = "deferred_exact"
    )
}

#' Build a validated pending custom metabolomics import generation
#'
#' @param stage Prepared import stage.
#' @param workflow_payload Established hydrated workflow payload.
#' @param manifest Exact assay manifest.
#' @param refs All immutable import refs.
#'
#' @return Validated payload-free pending generation.
#' @noRd
newMetabImportPendingStage <- function(stage, workflow_payload, manifest, refs) {
    spec <- metabArtifactImportSpec("custom")
    parameters <- list(
        capability_id = spec$capability_id,
        assay_order = manifest$assay_name,
        assay_roles = manifest$artifact_role,
        assay_identities = manifest$assay_identity,
        column_mapping = workflow_payload$columnMapping,
        column_mapping_serialized = artifactWorkflowStateSerializeMetadata(
            workflow_payload$columnMapping,
            "metabolomics import column mapping"
        ),
        readthrough_contract_version = 1L,
        input_format = "custom",
        data_level = "metabolite",
        sample_names_sanitized = isTRUE(
            workflow_payload$sampleNamesSanitized
        ),
        evidence_boundary = paste(
            "generated inputs certify scale only;",
            "reviewed fixtures certify parser and scientific parity"
        )
    )
    stage$refs <- refs
    pending <- list(
        schema = .METAB_IMPORT_PENDING_STAGE_SCHEMA,
        schema_version = 1L,
        stage = stage,
        parameters = parameters,
        source = metabArtifactAggregateSource(manifest),
        assay_manifest = manifest,
        process_evidence = list(
            writer_pid = as.integer(Sys.getpid()),
            verifier_pid = as.integer(Sys.getpid()),
            sequential_bounded_inline = TRUE,
            complete_source_materialized = FALSE,
            complete_encoded_payload_returned = FALSE
        )
    )
    validateMetabImportPendingStage(pending)
}

#' Validate one pending custom metabolomics import generation
#'
#' @param pending Candidate pending generation.
#'
#' @return The validated payload-free generation.
#' @noRd
validateMetabImportPendingStage <- function(pending) {
    valid <- is.list(pending) &&
        identical(pending$schema, .METAB_IMPORT_PENDING_STAGE_SCHEMA) &&
        identical(as.integer(pending$schema_version), 1L) &&
        is.list(pending$stage) && identical(pending$stage$stage_id, "import") &&
        is.list(pending$stage$refs) && length(pending$stage$refs) >= 4L &&
        isTRUE(pending$process_evidence$sequential_bounded_inline) &&
        identical(
            pending$process_evidence$complete_source_materialized,
            FALSE
        ) && identical(
            pending$process_evidence$complete_encoded_payload_returned,
            FALSE
        )
    if (!isTRUE(valid)) {
        metabArtifactAbort(
            "pending metabolomics import stage is malformed or unverified",
            "multischolar_invalid_metabolomics_pending_stage"
        )
    }
    artifactResourceDataOnly(pending, "metabolomics pending import stage")
    pending
}

#' Discard one unpublished custom metabolomics generation
#'
#' @param stage Pending stage bundle.
#'
#' @return `TRUE`, invisibly.
#' @noRd
discardMetabImportPendingStage <- function(stage) {
    store <- validateArtifactStore(stage$store)
    refs <- stage$refs %||% list()
    for (ref in refs) {
        normalized <- artifactStoreNormalizeRef(ref)
        managed <- artifactStoreManagedPaths(
            store,
            normalized$logical_key,
            normalized$artifact_id
        )
        payload_directory <- artifactStoreResolveFile(
            store,
            managed$payload_directory
        )
        sidecar <- artifactStoreResolveFile(store, managed$sidecar)
        intent <- artifactStoreResolveFile(store, managed$intent)
        unlink(payload_directory, recursive = TRUE, force = TRUE)
        unlink(c(sidecar, intent), force = TRUE)
    }
    invisible(TRUE)
}

#' Stage one exact custom metabolomics import
#'
#' @param workflow_data Mutable metabolomics workflow state.
#' @param assay1_file First assay source.
#' @param assay1_name First assay name.
#' @param assay2_file Optional second assay source.
#' @param assay2_name Optional second assay name.
#' @param column_mapping Established input mapping.
#' @param sanitize_names Whether sample names are sanitized.
#' @param resource_policy Optional project registry policy.
#'
#' @return Staged workflow payload and payload-free pending generation.
#' @noRd
stageMetabImportArtifacts <- function(
    workflow_data,
    assay1_file,
    assay1_name,
    assay2_file = NULL,
    assay2_name = NULL,
    column_mapping,
    sanitize_names = FALSE,
    resource_policy = NULL
) {
    prepared <- prepareMetabArtifactContext(
        workflow_data,
        format = "custom",
        data_type = "metabolite"
    )
    if (!isTRUE(prepared$enabled)) {
        return(list(
            enabled = FALSE,
            attempted = FALSE,
            ok = TRUE,
            reason = prepared$reason
        ))
    }
    sources <- metabImportStageSources(
        assay1_file,
        assay1_name,
        assay2_file,
        assay2_name
    )
    first_preflight <- metabImportStreamingPreflight(
        sources[[1L]],
        column_mapping,
        sanitize_names
    )
    output_mapping <- first_preflight$column_mapping
    stage <- metabImportPreparedStage(prepared, resource_policy)
    caller_owns_stage <- FALSE
    on.exit({
        if (!caller_owns_stage) {
            try(discardMetabImportPendingStage(stage), silent = TRUE)
        }
    }, add = TRUE)
    roles <- metabArtifactAssayRoles(names(sources))
    assay_refs <- lapply(seq_along(sources), \(index) {
        metabImportStageAssayRef(
            stage,
            roles[[index]],
            sources[[index]],
            column_mapping,
            sanitize_names
        )
    })
    names(assay_refs) <- unname(roles)
    proofs <- lapply(assay_refs, \(ref) {
        artifactStoreVerifyStreamingRef(
            stage$store,
            ref,
            verify_semantic = FALSE
        )
    })
    assays <- lapply(assay_refs, \(ref) {
        artifactStageReadTable(
            metabArtifactReadthroughAdapter(),
            stage$store,
            ref
        )
    })
    names(assays) <- names(sources)
    workflow_payload <- list(
        assayList = assays,
        sampleCols = output_mapping$sample_columns,
        sampleNamesSanitized = isTRUE(sanitize_names),
        dataFormat = "custom",
        columnMapping = output_mapping,
        processingLog = list(
            timestamp = Sys.time(),
            n_assays = length(assays),
            assay_names = names(assays),
            detected_format = "custom",
            n_metabolites = vapply(assays, \(assay) {
                metabolite_column <- output_mapping$metabolite_id_col
                if (metabolite_column %in% names(assay)) {
                    length(unique(assay[[metabolite_column]]))
                } else {
                    nrow(assay)
                }
            }, integer(1)),
            n_samples = length(output_mapping$sample_columns)
        ),
        sourceFiles = list(
            assay1 = unname(sources[[1L]]),
            assay2 = if (length(sources) > 1L) unname(sources[[2L]]) else NULL
        ),
        assaySourceFiles = as.list(sources)
    )
    spec <- metabArtifactImportSpec("custom")
    manifest <- metabArtifactImportManifest(
        workflow_payload,
        spec,
        sources,
        prepared$context$getIdentity()
    )
    metadata_tables <- metabArtifactImportTables(workflow_payload, manifest)[c(
        "assay_manifest", "column_mapping", "source_manifest"
    )]
    metadata_refs <- artifactStageWriteRefs(
        stage,
        metadata_tables,
        verification = "inline_exact",
        abort_fn = metabArtifactAbort
    )
    refs <- c(assay_refs, metadata_refs)
    pending <- newMetabImportPendingStage(
        stage,
        workflow_payload,
        manifest,
        refs
    )
    pending$process_evidence$proofs <- proofs
    result <- list(
        enabled = TRUE,
        attempted = TRUE,
        ok = TRUE,
        result = workflow_payload,
        pending_stage = pending,
        process_evidence = pending$process_evidence
    )
    caller_owns_stage <- TRUE
    result
}

#' Safely stage one exact custom metabolomics import
#'
#' @inheritParams stageMetabImportArtifacts
#' @param log_warn Warning logger used for memory fallback.
#'
#' @return Staged result or explicit memory fallback evidence.
#' @noRd
stageMetabImportArtifactsSafely <- function(
    workflow_data,
    assay1_file,
    assay1_name,
    assay2_file = NULL,
    assay2_name = NULL,
    column_mapping,
    sanitize_names = FALSE,
    resource_policy = NULL,
    log_warn = logger::log_warn
) {
    tryCatch(
        stageMetabImportArtifacts(
            workflow_data,
            assay1_file,
            assay1_name,
            assay2_file,
            assay2_name,
            column_mapping,
            sanitize_names,
            resource_policy
        ),
        error = \(error) {
            log_warn(paste(
                "metabolomics artifact import staging failed; using memory:",
                conditionMessage(error)
            ))
            list(
                enabled = TRUE,
                attempted = TRUE,
                ok = FALSE,
                error_class = class(error),
                error_message = conditionMessage(error)
            )
        }
    )
}

#' Publish one verified custom metabolomics import generation
#'
#' @param workflow_data Mutable metabolomics workflow state.
#' @param pending_stage Verified pending generation.
#'
#' @return Committed import stage evidence.
#' @noRd
publishMetabImportPendingStage <- function(workflow_data, pending_stage) {
    pending <- validateMetabImportPendingStage(pending_stage)
    context <- workflow_data$workflow_context
    if (!inherits(context, "WorkflowContext") || !context$isBound() ||
        !identical(context$getIdentity(), pending$stage$identity)) {
        metabArtifactAbort(
            "pending metabolomics import belongs to another workflow context",
            "multischolar_invalid_metabolomics_pending_context"
        )
    }
    registry <- projectRegistryForContext(
        context,
        pending$stage$resource_policy
    )
    session <- initializeProjectRegistry(registry)
    on.exit(closeProjectRegistry(session), add = TRUE)
    artifactWorkflowStateEnsureWorkflow(
        session,
        pending$stage$registry_identity
    )
    stage <- pending$stage
    stage$session <- session
    artifactStageRegister(
        stage,
        stage$refs,
        pending$parameters,
        pending$source,
        deferred_commit = FALSE
    )
    assay_refs <- stage$refs[pending$assay_manifest$artifact_role]
    names(assay_refs) <- pending$assay_manifest$assay_name
    list(
        enabled = TRUE,
        ok = TRUE,
        stage_id = "import",
        run_id = stage$run_id,
        action_id = stage$action_id,
        generation_id = stage$generation_id,
        refs = stage$refs,
        assay_refs = assay_refs,
        assay_manifest = pending$assay_manifest,
        committed = TRUE,
        process_evidence = pending$process_evidence
    )
}
