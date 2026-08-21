# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

protDiaEnrichRequestIdentifiers <- function(values) {
    values <- as.character(values)
    unique(values[!is.na(values) & nzchar(values)])
}

protDiaEnrichRequestBackground <- function(values, backend) {
    values <- as.character(values)
    if (identical(backend, "gprofiler2")) values <- values[!is.na(values)]
    unique(values)
}

protDiaEnrichIdentifierStats <- function(values, sent) {
    values <- as.character(values)
    valid <- values[!is.na(values) & nzchar(values)]
    list(
        input_count = as.integer(length(values)),
        missing_count = as.integer(sum(is.na(values))),
        blank_count = as.integer(sum(!is.na(values) & !nzchar(values))),
        duplicate_count = as.integer(length(valid) - length(unique(valid))),
        sent_count = as.integer(length(sent))
    )
}

protDiaEnrichCanonicalRequest <- function(
    backend,
    contrast,
    direction,
    identifiers,
    background,
    parameters,
    mapping = list()
) {
    backend <- match.arg(backend, c("gprofiler2", "clusterprofiler"))
    direction <- match.arg(direction, c("up", "down"))
    if (!is.list(mapping)) {
        protDiaEnrichArtifactAbort(
            "enrichment request mapping metadata must be a list",
            "multischolar_invalid_prot_dia_enrichment_request"
        )
    }
    sent_identifiers <- protDiaEnrichRequestIdentifiers(identifiers)
    sent_background <- protDiaEnrichRequestBackground(background, backend)
    mapping$identifier_stats <- protDiaEnrichIdentifierStats(
        identifiers,
        sent_identifiers
    )
    mapping$background_stats <- protDiaEnrichIdentifierStats(
        background,
        sent_background
    )
    request <- list(
        backend = backend,
        contrast = as.character(contrast),
        direction = direction,
        identifiers = sent_identifiers,
        background = sent_background,
        parameters = parameters,
        mapping = mapping
    )
    protDiaEnrichAssertSafeMetadata(
        request,
        "DIA-NN enrichment service request"
    )
    request
}

protDiaEnrichGprofilerContext <- function(
    contrast,
    direction,
    protein_id_column
) {
    list(
        contrast = contrast,
        direction = direction,
        mapping = list(
            protein_id_column = protein_id_column,
            transform = "remove_colon_suffix"
        )
    )
}

protDiaEnrichRecordGprofilerSkip <- function(
    execution_context,
    context,
    background,
    species,
    threshold,
    correction_method,
    exclude_iea
) {
    if (!isProtDiaEnrichExecutionContext(execution_context)) return(invisible(NULL))
    request <- protDiaEnrichCanonicalRequest(
        backend = "gprofiler2",
        contrast = context$contrast,
        direction = context$direction,
        identifiers = character(),
        background = background,
        parameters = list(
            organism = species,
            ordered_query = FALSE,
            sources = as.list(c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC")),
            user_threshold = as.double(threshold),
            correction_method = correction_method,
            exclude_iea = isTRUE(exclude_iea),
            evcodes = TRUE,
            domain_scope = "custom",
            significant = TRUE,
            highlight = TRUE
        ),
        mapping = context$mapping
    )
    protDiaEnrichRecordSkippedRequest(execution_context, request)
}

protDiaEnrichExecuteClusterRequest <- function(
    execution_context,
    genes,
    background,
    contrast,
    direction,
    taxon_id,
    q_cutoff,
    exclude_iea,
    protein_id_column,
    term2gene,
    term2name,
    enricher_fn
) {
    request <- protDiaEnrichCanonicalRequest(
        backend = "clusterprofiler",
        contrast = contrast,
        direction = direction,
        identifiers = genes,
        background = background,
        parameters = list(
            organism_taxid = as.character(taxon_id),
            sources = list("GO:BP", "GO:CC", "GO:MF"),
            user_threshold = as.double(q_cutoff),
            correction_method = "BH",
            method = "enricher",
            exclude_iea = isTRUE(exclude_iea)
        ),
        mapping = list(
            protein_id_column = protein_id_column,
            term2gene_digest = artifactSemanticDigest(term2gene),
            term2name_digest = artifactSemanticDigest(term2name),
            term2gene_rows = as.integer(nrow(term2gene)),
            term_count = as.integer(length(unique(term2gene$TERM))),
            mapped_identifier_count = as.integer(length(intersect(
                protDiaEnrichRequestIdentifiers(genes),
                unique(as.character(term2gene[[2L]]))
            ))),
            missing_mapping_count = as.integer(length(setdiff(
                protDiaEnrichRequestIdentifiers(genes),
                unique(as.character(term2gene[[2L]]))
            )))
        )
    )
    if (length(request$identifiers) == 0L) {
        protDiaEnrichRecordSkippedRequest(execution_context, request)
        return(NULL)
    }
    protDiaEnrichExecuteRequest(
        execution_context,
        request,
        call = function() enricher_fn(
            gene = request$identifiers,
            universe = request$background,
            TERM2GENE = term2gene,
            TERM2NAME = term2name,
            pvalueCutoff = q_cutoff,
            pAdjustMethod = "BH"
        )
    )
}

protDiaEnrichResponseTable <- function(backend, response) {
    if (identical(backend, "gprofiler2")) {
        if (!is.list(response) || !is.data.frame(response$result)) {
            protDiaEnrichArtifactAbort(
                "gprofiler2 returned an unsupported response",
                "multischolar_invalid_prot_dia_enrichment_response"
            )
        }
        return(response$result)
    }
    if (!methods::is(response, "enrichResult") &&
        !(isS4(response) && "result" %in% methods::slotNames(response))) {
        protDiaEnrichArtifactAbort(
            "clusterProfiler returned an unsupported response",
            "multischolar_invalid_prot_dia_enrichment_response"
        )
    }
    response@result
}

protDiaEnrichSafeGprofilerMeta <- function(response, request) {
    query_metadata <- tryCatch(
        response$meta$query_metadata,
        error = function(error) NULL
    )
    threshold <- query_metadata$user_threshold
    if (is.null(threshold)) threshold <- request$parameters$user_threshold
    list(
        query_metadata = list(
            user_threshold = as.double(threshold),
            organism = as.character(request$parameters$organism),
            sources = as.list(request$parameters$sources),
            correction_method = as.character(
                request$parameters$correction_method
            ),
            domain_scope = as.character(request$parameters$domain_scope),
            exclude_iea = isTRUE(request$parameters$exclude_iea)
        )
    )
}

protDiaEnrichSafeResponseMeta <- function(response, request) {
    if (identical(request$backend, "gprofiler2")) {
        return(protDiaEnrichSafeGprofilerMeta(response, request))
    }
    list(
        local_execution = list(
            method = as.character(request$parameters$method),
            organism_taxid = as.character(request$parameters$organism_taxid),
            sources = as.list(request$parameters$sources),
            user_threshold = as.double(request$parameters$user_threshold),
            correction_method = as.character(
                request$parameters$correction_method
            ),
            exclude_iea = isTRUE(request$parameters$exclude_iea)
        )
    )
}

protDiaEnrichValidateIdentifierStats <- function(value, owner) {
    required <- c(
        "input_count", "missing_count", "blank_count", "duplicate_count",
        "sent_count"
    )
    valid <- is.list(value) && identical(names(value), required) &&
        all(vapply(value, function(item) {
            length(item) == 1L && is.numeric(item) && !is.na(item) &&
                is.finite(item) && item >= 0 && item == floor(item)
        }, logical(1)))
    if (!isTRUE(valid)) {
        protDiaEnrichArtifactAbort(
            sprintf("%s identifier statistics are malformed", owner),
            "multischolar_invalid_prot_dia_enrichment_request"
        )
    }
    lapply(value, as.integer)
}

protDiaEnrichValidateRequestParameters <- function(value, backend) {
    required <- if (identical(backend, "gprofiler2")) {
        c(
            "organism", "ordered_query", "sources", "user_threshold",
            "correction_method", "exclude_iea", "evcodes", "domain_scope",
            "significant", "highlight"
        )
    } else {
        c(
            "organism_taxid", "sources", "user_threshold",
            "correction_method", "method", "exclude_iea"
        )
    }
    valid <- is.list(value) && identical(names(value), required) &&
        is.list(value$sources) && length(value$sources) > 0L &&
        workflowCapabilityScalarString(value$correction_method) &&
        length(value$user_threshold) == 1L &&
        is.numeric(value$user_threshold) && is.finite(value$user_threshold)
    if (identical(backend, "gprofiler2")) {
        valid <- valid && workflowCapabilityScalarString(value$organism) &&
            workflowCapabilityScalarString(value$domain_scope) &&
            all(vapply(
                value[c(
                    "ordered_query", "exclude_iea", "evcodes", "significant",
                    "highlight"
                )],
                function(item) length(item) == 1L && is.logical(item),
                logical(1)
            ))
    } else {
        valid <- valid && workflowCapabilityScalarString(value$organism_taxid) &&
            identical(value$method, "enricher") &&
            identical(value$correction_method, "BH") &&
            length(value$exclude_iea) == 1L && is.logical(value$exclude_iea)
    }
    if (!isTRUE(valid)) {
        protDiaEnrichArtifactAbort(
            "enrichment request parameters are malformed",
            "multischolar_invalid_prot_dia_enrichment_request"
        )
    }
    value$user_threshold <- as.double(value$user_threshold)
    value
}

protDiaEnrichValidateRequestMapping <- function(value, backend) {
    shared <- c("identifier_stats", "background_stats")
    required <- if (identical(backend, "gprofiler2")) {
        c("protein_id_column", "transform", shared)
    } else {
        c(
            "protein_id_column", "term2gene_digest", "term2name_digest",
            "term2gene_rows", "term_count", "mapped_identifier_count",
            "missing_mapping_count", shared
        )
    }
    valid <- is.list(value) && identical(names(value), required) &&
        workflowCapabilityScalarString(value$protein_id_column)
    if (identical(backend, "gprofiler2")) {
        valid <- valid && identical(value$transform, "remove_colon_suffix")
    } else {
        counts <- value[c(
            "term2gene_rows", "term_count", "mapped_identifier_count",
            "missing_mapping_count"
        )]
        valid <- valid && all(vapply(counts, function(item) {
            length(item) == 1L && is.numeric(item) && !is.na(item) &&
                is.finite(item) && item >= 0 && item == floor(item)
        }, logical(1)))
        artifactRefValidateDigest(value$term2gene_digest, "TERM2GENE digest")
        artifactRefValidateDigest(value$term2name_digest, "TERM2NAME digest")
        value[names(counts)] <- lapply(counts, as.integer)
    }
    if (!isTRUE(valid)) {
        protDiaEnrichArtifactAbort(
            "enrichment request mapping metadata is malformed",
            "multischolar_invalid_prot_dia_enrichment_request"
        )
    }
    value$identifier_stats <- protDiaEnrichValidateIdentifierStats(
        value$identifier_stats,
        "query"
    )
    value$background_stats <- protDiaEnrichValidateIdentifierStats(
        value$background_stats,
        "background"
    )
    value
}

protDiaEnrichValidateRequestMetadata <- function(value, backend) {
    required <- c(
        "backend", "contrast", "direction", "query_count", "query_digest",
        "background_count", "background_digest", "parameters", "mapping"
    )
    valid <- is.list(value) && identical(names(value), required) &&
        identical(value$backend, backend) &&
        workflowCapabilityScalarString(value$contrast) &&
        value$direction %in% c("up", "down")
    counts <- value[c("query_count", "background_count")]
    valid <- valid && all(vapply(counts, function(item) {
        length(item) == 1L && is.numeric(item) && !is.na(item) &&
            is.finite(item) && item >= 0 && item == floor(item)
    }, logical(1)))
    if (!isTRUE(valid)) {
        protDiaEnrichArtifactAbort(
            "enrichment request metadata is malformed",
            "multischolar_invalid_prot_dia_enrichment_request"
        )
    }
    artifactRefValidateDigest(value$query_digest, "enrichment query digest")
    artifactRefValidateDigest(
        value$background_digest,
        "enrichment background digest"
    )
    value$query_count <- as.integer(value$query_count)
    value$background_count <- as.integer(value$background_count)
    value$parameters <- protDiaEnrichValidateRequestParameters(
        value$parameters,
        backend
    )
    value$mapping <- protDiaEnrichValidateRequestMapping(value$mapping, backend)
    value
}

protDiaEnrichValidateSoftware <- function(value, backend) {
    required <- c(
        "multischolar", "backend", "backend_package", "backend_version", "r"
    )
    expected_package <- if (identical(backend, "gprofiler2")) {
        "gprofiler2"
    } else {
        "clusterProfiler"
    }
    valid <- is.list(value) && identical(names(value), required) &&
        all(vapply(value, workflowCapabilityScalarString, logical(1))) &&
        identical(value$backend, backend) &&
        identical(value$backend_package, expected_package)
    if (!isTRUE(valid)) {
        protDiaEnrichArtifactAbort(
            "enrichment software provenance is malformed",
            "multischolar_invalid_prot_dia_enrichment_manifest"
        )
    }
    value
}

protDiaEnrichValidateResponseMeta <- function(value, request) {
    if (identical(request$backend, "gprofiler2")) {
        valid <- is.list(value) && identical(names(value), "query_metadata") &&
            is.list(value$query_metadata)
        metadata <- value$query_metadata
        required <- c(
            "user_threshold", "organism", "sources", "correction_method",
            "domain_scope", "exclude_iea"
        )
        valid <- valid && identical(names(metadata), required) &&
            identical(as.double(metadata$user_threshold),
                      request$parameters$user_threshold) &&
            identical(metadata$organism, request$parameters$organism) &&
            identical(metadata$sources, request$parameters$sources) &&
            identical(
                metadata$correction_method,
                request$parameters$correction_method
            ) && identical(
                metadata$domain_scope,
                request$parameters$domain_scope
            ) && identical(
                metadata$exclude_iea,
                request$parameters$exclude_iea
            )
    } else {
        valid <- is.list(value) && identical(names(value), "local_execution") &&
            is.list(value$local_execution)
        metadata <- value$local_execution
        required <- c(
            "method", "organism_taxid", "sources", "user_threshold",
            "correction_method", "exclude_iea"
        )
        valid <- valid && identical(names(metadata), required) &&
            identical(metadata$method, request$parameters$method) &&
            identical(
                metadata$organism_taxid,
                request$parameters$organism_taxid
            ) && identical(metadata$sources, request$parameters$sources) &&
            identical(as.double(metadata$user_threshold),
                      request$parameters$user_threshold) &&
            identical(
                metadata$correction_method,
                request$parameters$correction_method
            ) && identical(
                metadata$exclude_iea,
                request$parameters$exclude_iea
            )
    }
    if (!isTRUE(valid)) {
        protDiaEnrichArtifactAbort(
            "enrichment response metadata differs from its request",
            "multischolar_prot_dia_enrichment_request_mismatch"
        )
    }
    value
}

protDiaEnrichRequestFromManifest <- function(metadata, identifiers) {
    query_rows <- which(identifiers$role == "query")
    background_rows <- which(identifiers$role == "background")
    valid <- all(identifiers$role %in% c("query", "background")) &&
        identical(identifiers$sequence[query_rows], seq_along(query_rows)) &&
        identical(
            identifiers$sequence[background_rows],
            seq_along(background_rows)
        )
    if (!isTRUE(valid)) {
        protDiaEnrichArtifactAbort(
            "enrichment identifier order is malformed",
            "multischolar_prot_dia_enrichment_request_mismatch"
        )
    }
    request <- list(
        backend = metadata$backend,
        contrast = metadata$contrast,
        direction = metadata$direction,
        identifiers = identifiers$identifier[query_rows],
        background = identifiers$identifier[background_rows],
        parameters = metadata$parameters,
        mapping = metadata$mapping
    )
    valid <- identical(length(request$identifiers), metadata$query_count) &&
        identical(length(request$background), metadata$background_count) &&
        identical(
            artifactSemanticDigest(request$identifiers),
            metadata$query_digest
        ) && identical(
            artifactSemanticDigest(request$background),
            metadata$background_digest
        )
    if (!isTRUE(valid)) {
        protDiaEnrichArtifactAbort(
            "enrichment identifier table differs from its request metadata",
            "multischolar_prot_dia_enrichment_request_mismatch"
        )
    }
    request
}

protDiaEnrichBuildResponse <- function(manifest, table) {
    if (identical(manifest$backend, "gprofiler2")) {
        return(list(result = table, meta = manifest$response_meta))
    }
    protDiaEnrichArtifactAbort(
        "clusterProfiler responses are local-only and cannot be replayed",
        "multischolar_unsupported_prot_dia_enrichment_replay"
    )
}

protDiaEnrichExecuteWithRetries <- function(
    call,
    max_retries,
    wait_time,
    sleep_fn = Sys.sleep
) {
    last_error <- NULL
    for (attempt in seq_len(as.integer(max_retries))) {
        result <- tryCatch(
            call(),
            error = function(error) {
                last_error <<- error
                NULL
            }
        )
        if (!is.null(result)) {
            return(list(result = result, attempts = attempt, error = NULL))
        }
        if (attempt < max_retries && wait_time > 0) sleep_fn(wait_time)
    }
    list(result = NULL, attempts = as.integer(max_retries), error = last_error)
}

newProtDiaEnrichExecutionContext <- function(
    workflow_data,
    mode = getOption("multischolar.prot_dia.enrichment_service_mode", "auto"),
    sleep_fn = Sys.sleep,
    failure_injector = NULL
) {
    mode <- match.arg(mode, c("auto", "live", "replay"))
    records <- list()
    cache <- new.env(parent = emptyenv())

    append_record <- function(record) {
        records[[length(records) + 1L]] <<- record
        invisible(record)
    }

    append_failed <- function(request, request_id, request_digest, state, attempts) {
        append_record(list(
            request_id = request_id,
            request_digest = request_digest,
            backend = request$backend,
            contrast = request$contrast,
            direction = request$direction,
            status = "failed",
            execution_state = state,
            attempts = as.integer(attempts),
            response = NULL
        ))
    }

    execute <- function(
        request,
        call,
        max_retries = 1L,
        wait_time = 0
    ) {
        request_digest <- artifactSemanticDigest(request)
        request_id <- protDiaEnrichRequestId(request)
        cached <- if (identical(request$backend, "gprofiler2")) {
            cache[[request_id]]
        } else {
            NULL
        }
        if (!is.null(cached)) {
            record <- cached$record
            record$execution_state <- "cache"
            append_record(record)
            return(cached$result)
        }

        if (identical(request$backend, "gprofiler2") &&
            !identical(mode, "live")) {
            replay <- protDiaEnrichLoadServiceResponse(
                workflow_data,
                request,
                required = identical(mode, "replay")
            )
            if (!is.null(replay)) {
                record <- replay$record
                record$execution_state <- "replay"
                append_record(record)
                cache[[request_id]] <- list(
                    result = replay$result,
                    record = record
                )
                return(replay$result)
            }
        }

        if (identical(mode, "replay") &&
            identical(request$backend, "gprofiler2")) {
            protDiaEnrichArtifactAbort(
                "required enrichment response replay is unavailable",
                "multischolar_missing_prot_dia_enrichment_replay"
            )
        }

        outcome <- protDiaEnrichExecuteWithRetries(
            call,
            max_retries,
            wait_time,
            sleep_fn
        )
        state <- if (identical(request$backend, "gprofiler2")) {
            "live"
        } else {
            "local"
        }
        if (is.null(outcome$result)) {
            append_failed(
                request,
                request_id,
                request_digest,
                state,
                outcome$attempts
            )
            stop(if (is.null(outcome$error)) {
                "enrichment service returned NULL after all attempts"
            } else {
                conditionMessage(outcome$error)
            }, call. = FALSE)
        }
        response <- tryCatch(
            protDiaEnrichPersistServiceResponse(
                workflow_data,
                request,
                outcome$result,
                failure_injector
            ),
            error = function(error) {
                append_failed(
                    request,
                    request_id,
                    request_digest,
                    state,
                    outcome$attempts
                )
                stop(error)
            }
        )
        record <- list(
            request_id = request_id,
            request_digest = request_digest,
            backend = request$backend,
            contrast = request$contrast,
            direction = request$direction,
            status = "succeeded",
            execution_state = state,
            attempts = outcome$attempts,
            response = response[names(response) != "result"]
        )
        append_record(record)
        if (identical(request$backend, "gprofiler2")) {
            cache[[request_id]] <- list(result = outcome$result, record = record)
        }
        outcome$result
    }

    record_skip <- function(request) {
        append_record(list(
            request_id = protDiaEnrichRequestId(request),
            request_digest = artifactSemanticDigest(request),
            backend = request$backend,
            contrast = request$contrast,
            direction = request$direction,
            status = "skipped_empty_input",
            execution_state = "not_called",
            attempts = 0L,
            response = NULL
        ))
        invisible(NULL)
    }

    structure(
        list(
            execute = execute,
            recordSkip = record_skip,
            getRecords = function() records,
            assertComplete = function() {
                failed <- vapply(
                    records,
                    function(record) identical(record$status, "failed"),
                    logical(1)
                )
                if (any(failed)) {
                    protDiaEnrichArtifactAbort(
                        "one or more enrichment service calls failed",
                        "multischolar_prot_dia_enrichment_service_failed"
                    )
                }
                invisible(TRUE)
            }
        ),
        class = c("MultiScholaRProtDiaEnrichExecution", "list")
    )
}

isProtDiaEnrichExecutionContext <- function(value) {
    inherits(value, "MultiScholaRProtDiaEnrichExecution") &&
        is.function(value$execute) && is.function(value$getRecords)
}

protDiaEnrichExecuteRequest <- function(
    execution_context,
    request,
    call,
    max_retries = 1L,
    wait_time = 0
) {
    if (!isProtDiaEnrichExecutionContext(execution_context)) return(call())
    execution_context$execute(request, call, max_retries, wait_time)
}

protDiaEnrichRecordSkippedRequest <- function(execution_context, request) {
    if (isProtDiaEnrichExecutionContext(execution_context)) {
        execution_context$recordSkip(request)
    }
    invisible(NULL)
}

protDiaEnrichValidateExecutionRecords <- function(
    execution_context,
    backend,
    contrast
) {
    if (!isProtDiaEnrichExecutionContext(execution_context)) {
        protDiaEnrichArtifactAbort(
            "artifact enrichment requires its execution provenance context",
            "multischolar_missing_prot_dia_enrichment_provenance"
        )
    }
    execution_context$assertComplete()
    records <- execution_context$getRecords()
    valid <- length(records) == 2L &&
        setequal(vapply(records, `[[`, character(1), "direction"), c("up", "down")) &&
        all(vapply(records, function(record) {
            identical(record$backend, backend) &&
                identical(record$contrast, contrast)
        }, logical(1)))
    if (!isTRUE(valid)) {
        protDiaEnrichArtifactAbort(
            "enrichment provenance is incomplete for the selected contrast",
            "multischolar_missing_prot_dia_enrichment_provenance"
        )
    }
    records
}
