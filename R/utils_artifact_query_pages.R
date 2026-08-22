# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

.ARTIFACT_QUERY_PAGE_SCHEMA <- "multischolar.artifact_query_page_specification"
.ARTIFACT_QUERY_PAGE_VERSION <- 1L
.ARTIFACT_QUERY_CURSOR_SCHEMA <- "multischolar.artifact_query_cursor"
.ARTIFACT_QUERY_CURSOR_VERSION <- 1L

artifactQueryPageScalarInteger <- function(value, field, maximum = 2^31 - 1L) {
    valid <- length(value) == 1L && is.numeric(value) && !is.na(value) &&
        is.finite(value) && value >= 1 && value == floor(value) &&
        value <= maximum
    if (!isTRUE(valid)) {
        artifactQueryAbort(
            sprintf("artifact page field '%s' must be a bounded integer", field),
            "multischolar_invalid_artifact_query_page"
        )
    }
    as.integer(value)
}

validateArtifactQueryPageSpecification <- function(specification) {
    required <- c(
        "schema", "schema_version", "query_id", "state_role",
        "projections", "filters", "sorts", "default_sort",
        "max_rows", "max_bytes"
    )
    valid <- is.list(specification) && identical(names(specification), required) &&
        identical(specification$schema, .ARTIFACT_QUERY_PAGE_SCHEMA) &&
        identical(
            workflowStateVersionValue(specification$schema_version),
            .ARTIFACT_QUERY_PAGE_VERSION
        ) &&
        workflowCapabilityScalarString(specification$query_id) &&
        workflowCapabilityScalarString(specification$state_role)
    if (!isTRUE(valid)) {
        artifactQueryAbort(
            "artifact page specification is malformed or unsupported",
            "multischolar_invalid_artifact_query_page"
        )
    }
    projections <- specification$projections
    if (is.list(projections) &&
        all(vapply(projections, is.character, logical(1)))) {
        projections <- unlist(projections, use.names = FALSE)
    }
    if (!is.character(projections) || length(projections) == 0L ||
        anyNA(projections) || any(!nzchar(projections)) ||
        anyDuplicated(projections) > 0L) {
        artifactQueryAbort(
            "artifact page projections must be unique logical columns",
            "multischolar_invalid_artifact_query_page"
        )
    }
    filters <- specification$filters
    if (!is.list(filters) || (length(filters) > 0L && is.null(names(filters))) ||
        any(!nzchar(names(filters))) || anyDuplicated(names(filters)) > 0L) {
        artifactQueryAbort(
            "artifact page filters must be an unambiguous named list",
            "multischolar_invalid_artifact_query_page"
        )
    }
    allowed_operators <- c(
        "equal", "in", "gt", "gte", "lt", "lte", "between",
        "is_missing", "contains"
    )
    for (filter_id in names(filters)) {
        filter <- filters[[filter_id]]
        if (is.list(filter$operators)) {
            filter$operators <- unlist(filter$operators, use.names = FALSE)
        }
        filter_valid <- is.list(filter) &&
            identical(names(filter), c("column", "type", "operators")) &&
            workflowCapabilityScalarString(filter$column) &&
            filter$column %in% projections &&
            workflowCapabilityScalarString(filter$type) &&
            filter$type %in% c("character", "integer", "double", "logical") &&
            is.character(filter$operators) && length(filter$operators) > 0L &&
            !anyNA(filter$operators) &&
            all(filter$operators %in% allowed_operators) &&
            anyDuplicated(filter$operators) == 0L
        if (!isTRUE(filter_valid)) {
            artifactQueryAbort(
                sprintf("artifact page filter '%s' is malformed", filter_id),
                "multischolar_invalid_artifact_query_page"
            )
        }
        if ("contains" %in% filter$operators &&
            !identical(filter$type, "character")) {
            artifactQueryAbort(
                "artifact page contains filters require character columns",
                "multischolar_invalid_artifact_query_page"
            )
        }
        filters[[filter_id]] <- filter
    }
    sorts <- specification$sorts
    if (!is.list(sorts) || length(sorts) == 0L || is.null(names(sorts)) ||
        any(!nzchar(names(sorts))) || anyDuplicated(names(sorts)) > 0L) {
        artifactQueryAbort(
            "artifact page sorts must be an unambiguous named list",
            "multischolar_invalid_artifact_query_page"
        )
    }
    for (sort_id in names(sorts)) {
        sort <- sorts[[sort_id]]
        if (is.list(sort$directions)) {
            sort$directions <- unlist(sort$directions, use.names = FALSE)
        }
        sort_valid <- is.list(sort) &&
            identical(names(sort), c("column", "transform", "directions")) &&
            workflowCapabilityScalarString(sort$column) &&
            sort$column %in% projections &&
            workflowCapabilityScalarString(sort$transform) &&
            sort$transform %in% c("identity", "absolute") &&
            is.character(sort$directions) && length(sort$directions) > 0L &&
            !anyNA(sort$directions) &&
            all(sort$directions %in% c("asc", "desc")) &&
            anyDuplicated(sort$directions) == 0L
        if (!isTRUE(sort_valid)) {
            artifactQueryAbort(
                sprintf("artifact page sort '%s' is malformed", sort_id),
                "multischolar_invalid_artifact_query_page"
            )
        }
        sorts[[sort_id]] <- sort
    }
    default_sort <- specification$default_sort
    if (!is.list(default_sort) ||
        !identical(names(default_sort), c("sort_id", "direction")) ||
        !workflowCapabilityScalarString(default_sort$sort_id) ||
        !default_sort$sort_id %in% names(sorts) ||
        !workflowCapabilityScalarString(default_sort$direction) ||
        !default_sort$direction %in% sorts[[default_sort$sort_id]]$directions) {
        artifactQueryAbort(
            "artifact page default sort is not declared",
            "multischolar_invalid_artifact_query_page"
        )
    }
    specification$schema_version <- .ARTIFACT_QUERY_PAGE_VERSION
    specification$projections <- projections
    specification$filters <- filters
    specification$sorts <- sorts
    specification$max_rows <- artifactQueryPageScalarInteger(
        specification$max_rows, "max_rows"
    )
    specification$max_bytes <- artifactQueryPageScalarInteger(
        specification$max_bytes, "max_bytes"
    )
    specification
}

newArtifactQueryPageSpecification <- function(
    query_id,
    state_role,
    projections,
    filters,
    sorts,
    default_sort,
    max_rows,
    max_bytes
) {
    validateArtifactQueryPageSpecification(list(
        schema = .ARTIFACT_QUERY_PAGE_SCHEMA,
        schema_version = .ARTIFACT_QUERY_PAGE_VERSION,
        query_id = query_id,
        state_role = state_role,
        projections = projections,
        filters = filters,
        sorts = sorts,
        default_sort = default_sort,
        max_rows = max_rows,
        max_bytes = max_bytes
    ))
}

artifactQueryPageReference <- function(store, ref, specification) {
    store <- validateArtifactStore(store)
    ref <- artifactStoreNormalizeRef(ref)
    owned <- identical(ref$logical_key$project_id, store$project_id) &&
        identical(ref$logical_key$omic_type, store$labels$omic_type) &&
        identical(
            ref$logical_key$workflow_slug,
            store$labels$workflow_slug
        ) &&
        identical(ref$logical_key$state_role, specification$state_role)
    if (!isTRUE(owned)) {
        artifactQueryAbort(
            "artifact page query does not own this project or state role",
            "multischolar_artifact_query_role_mismatch"
        )
    }
    sidecar_path <- artifactStoreManagedPaths(
        store, ref$logical_key, ref$artifact_id
    )$sidecar
    sidecar <- artifactStoreReadSidecar(
        store, sidecar_path, validate_payload = FALSE
    )
    if (!identical(sidecar$artifact_ref, ref)) {
        artifactQueryAbort(
            "artifact page reference differs from its immutable sidecar",
            "multischolar_artifact_query_ref_mismatch"
        )
    }
    payload_path <- artifactStoreResolveFile(
        store, ref$relative_path, must_exist = TRUE
    )
    bytes <- unname(as.numeric(file.info(payload_path)$size))
    if (!identical(bytes, as.numeric(ref$shape$bytes))) {
        artifactQueryAbort(
            "artifact page payload size differs from its immutable reference",
            "multischolar_artifact_query_payload_mismatch"
        )
    }
    list(
        ref = ref,
        sidecar = sidecar,
        payload_path = payload_path
    )
}

artifactQueryPageSort <- function(specification, sort_id, direction, columns) {
    if (is.null(sort_id)) sort_id <- specification$default_sort$sort_id
    if (is.null(direction)) direction <- specification$default_sort$direction
    if (!workflowCapabilityScalarString(sort_id) ||
        !sort_id %in% names(specification$sorts)) {
        artifactQueryAbort(
            "artifact page sort is not declared",
            "multischolar_invalid_artifact_query_sort"
        )
    }
    declaration <- specification$sorts[[sort_id]]
    if (!workflowCapabilityScalarString(direction) ||
        !direction %in% declaration$directions) {
        artifactQueryAbort(
            "artifact page sort direction is not declared",
            "multischolar_invalid_artifact_query_sort"
        )
    }
    descriptor <- columns$descriptors[[declaration$column]]
    if (is.null(descriptor)) {
        artifactQueryAbort(
            "artifact page sort column is absent from the payload",
            "multischolar_artifact_query_schema_mismatch"
        )
    }
    if (identical(declaration$transform, "absolute") &&
        !descriptor$r_type %in% c("integer", "integer64", "double")) {
        artifactQueryAbort(
            "artifact absolute sort requires a numeric payload column",
            "multischolar_artifact_query_schema_mismatch"
        )
    }
    list(sort_id = sort_id, direction = direction, declaration = declaration)
}

artifactQueryPageDefaultSql <- function(type) {
    if (type %in% c("integer", "integer64", "double", "Date", "POSIXct")) {
        return("0")
    }
    if (identical(type, "logical")) return("FALSE")
    "''"
}

artifactQueryPageOrderTerms <- function(connection, sort, columns) {
    quote <- function(value) {
        as.character(DBI::dbQuoteIdentifier(connection, value))
    }
    descriptor <- columns$descriptors[[sort$declaration$column]]
    physical <- quote(descriptor$physical_name)
    status <- if (is.null(descriptor$status_name)) {
        paste0("CASE WHEN ", physical, " IS NULL THEN 1 ELSE 0 END")
    } else {
        quote(descriptor$status_name)
    }
    value <- paste0(
        "COALESCE(", physical, ", ",
        artifactQueryPageDefaultSql(descriptor$r_type), ")"
    )
    if (identical(sort$declaration$transform, "absolute")) {
        value <- paste0("ABS(", value, ")")
    }
    list(
        list(expression = status, direction = "asc", type = "integer"),
        list(
            expression = value,
            direction = sort$direction,
            type = descriptor$r_type
        ),
        list(
            expression = quote(.artifactRowOrderColumn),
            direction = "asc",
            type = "double"
        )
    )
}

artifactQueryPageValidateCursor <- function(cursor, specification, ref, sort) {
    if (is.null(cursor)) return(NULL)
    required <- c(
        "schema", "schema_version", "query_id", "artifact_id",
        "sort_id", "direction", "values"
    )
    valid <- is.list(cursor) && identical(names(cursor), required) &&
        identical(cursor$schema, .ARTIFACT_QUERY_CURSOR_SCHEMA) &&
        identical(
            workflowStateVersionValue(cursor$schema_version),
            .ARTIFACT_QUERY_CURSOR_VERSION
        ) &&
        identical(cursor$query_id, specification$query_id) &&
        identical(cursor$artifact_id, ref$artifact_id) &&
        identical(cursor$sort_id, sort$sort_id) &&
        identical(cursor$direction, sort$direction) &&
        is.list(cursor$values) && length(cursor$values) == 3L &&
        all(vapply(cursor$values, function(value) {
            length(value) == 1L && !is.na(value) &&
                (is.character(value) || is.logical(value) ||
                    (is.numeric(value) && is.finite(value)))
        }, logical(1)))
    if (!isTRUE(valid)) {
        artifactQueryAbort(
            "artifact page cursor is malformed, stale, or operation-mismatched",
            "multischolar_invalid_artifact_query_cursor"
        )
    }
    cursor
}

artifactQueryPageCursorSql <- function(terms, cursor) {
    if (is.null(cursor)) return(list(sql = "", values = list()))
    clauses <- character()
    values <- list()
    for (index in seq_along(terms)) {
        prefix <- character()
        if (index > 1L) {
            prefix <- vapply(seq_len(index - 1L), function(previous) {
                paste0(terms[[previous]]$expression, " = ?")
            }, character(1))
            values <- c(values, cursor$values[seq_len(index - 1L)])
        }
        comparator <- if (identical(terms[[index]]$direction, "desc")) {
            " < ?"
        } else {
            " > ?"
        }
        clauses <- c(clauses, paste0(
            "(",
            paste(
                c(prefix, paste0(terms[[index]]$expression, comparator)),
                collapse = " AND "
            ),
            ")"
        ))
        values <- c(values, cursor$values[index])
    }
    list(sql = paste0("(", paste(clauses, collapse = " OR "), ")"), values = values)
}

artifactQueryPageStatements <- function(
    connection,
    projections,
    columns,
    filter_sql,
    cursor_sql,
    terms
) {
    quote <- function(value) {
        as.character(DBI::dbQuoteIdentifier(connection, value))
    }
    projection_terms <- paste0(
        quote(columns$physical[projections]), " AS ", quote(projections)
    )
    status_projections <- projections[vapply(
        columns$descriptors[projections],
        function(descriptor) !is.null(descriptor$status_name),
        logical(1)
    )]
    status_aliases <- stats::setNames(
        sprintf(".multischolar_page_status_%04d", seq_along(status_projections)),
        status_projections
    )
    status_terms <- vapply(status_projections, function(projection) {
        paste0(
            quote(columns$descriptors[[projection]]$status_name),
            " AS ", quote(status_aliases[[projection]])
        )
    }, character(1))
    cursor_aliases <- sprintf(".multischolar_cursor_%04d", seq_along(terms))
    cursor_terms <- vapply(seq_along(terms), function(index) {
        paste0(terms[[index]]$expression, " AS ", quote(cursor_aliases[[index]]))
    }, character(1))
    where <- sub("^ WHERE ", "", filter_sql)
    where <- where[nzchar(where)]
    if (nzchar(cursor_sql)) where <- c(where[nzchar(where)], cursor_sql)
    where <- if (length(where) == 0L) "" else {
        paste0(" WHERE ", paste(where, collapse = " AND "))
    }
    order <- paste(vapply(terms, function(term) {
        paste(term$expression, toupper(term$direction))
    }, character(1)), collapse = ", ")
    query <- paste0(
        "SELECT ",
        paste(c(projection_terms, status_terms, cursor_terms), collapse = ", "),
        " FROM read_parquet(?)", where, " ORDER BY ", order, " LIMIT ?"
    )
    byte_terms <- paste0(
        "octet_length(encode(COALESCE(CAST(", quote(projections),
        " AS VARCHAR), '')))"
    )
    estimate <- paste0(
        "WITH bounded_artifact_page AS (", query, ") ",
        "SELECT COUNT(*) AS row_count, COALESCE(SUM(",
        paste(byte_terms, collapse = " + "),
        "), 0) AS payload_bytes FROM bounded_artifact_page"
    )
    list(
        query = query,
        estimate = estimate,
        status_aliases = status_aliases,
        cursor_aliases = cursor_aliases
    )
}

queryArtifactPage <- function(
    store,
    ref,
    specification,
    projections = NULL,
    filters = list(),
    sort_id = NULL,
    direction = NULL,
    cursor = NULL,
    limit = NULL,
    resource_policy = NULL,
    query_session = NULL
) {
    specification <- validateArtifactQueryPageSpecification(specification)
    source <- artifactQueryPageReference(store, ref, specification)
    columns <- artifactQueryColumnMap(source$sidecar$codec_metadata)
    projections <- artifactQueryNormalizeProjections(
        specification, projections, columns
    )
    bounds <- artifactQueryNormalizeBounds(specification, limit, resource_policy)
    artifactQueryAssertRss(bounds$policy, "before_page_query")
    lease <- artifactQueryLease(store, bounds$policy, query_session)
    handle <- lease$handle
    if (isTRUE(lease$transient)) {
        on.exit(artifactQueryDisconnect(handle), add = TRUE)
    }
    sort <- artifactQueryPageSort(
        specification, sort_id, direction, columns
    )
    cursor <- artifactQueryPageValidateCursor(
        cursor, specification, source$ref, sort
    )
    filter <- artifactQueryFilterSql(
        specification, filters, columns, handle$connection
    )
    terms <- artifactQueryPageOrderTerms(handle$connection, sort, columns)
    cursor_filter <- artifactQueryPageCursorSql(terms, cursor)
    statements <- artifactQueryPageStatements(
        handle$connection,
        projections,
        columns,
        filter$sql,
        cursor_filter$sql,
        terms
    )
    fetch_limit <- as.integer(bounds$limit + 1L)
    values <- c(
        list(source$payload_path),
        filter$values,
        cursor_filter$values,
        list(fetch_limit)
    )
    estimate <- tryCatch(
        projectRegistryFetchBound(
            handle$connection, statements$estimate, values
        ),
        error = function(error) artifactQueryAbort(
            "bounded artifact page byte estimation failed",
            "multischolar_artifact_query_failed",
            parent = error
        )
    )
    estimated_bytes <- as.numeric(estimate$payload_bytes[[1L]]) * 2 +
        as.numeric(estimate$row_count[[1L]]) * 256
    if (estimated_bytes > bounds$max_bytes) {
        artifactQueryAbort(
            "artifact page exceeds its pre-materialization byte limit",
            "multischolar_artifact_query_byte_limit_exceeded",
            estimated_bytes = estimated_bytes
        )
    }
    result <- tryCatch(
        projectRegistryFetchBound(handle$connection, statements$query, values),
        error = function(error) artifactQueryAbort(
            "bounded artifact page query failed",
            "multischolar_artifact_query_failed",
            parent = error
        )
    )
    has_more <- nrow(result) > bounds$limit
    if (has_more) result <- result[seq_len(bounds$limit), , drop = FALSE]
    next_cursor <- NULL
    if (has_more && nrow(result) > 0L) {
        next_cursor <- list(
            schema = .ARTIFACT_QUERY_CURSOR_SCHEMA,
            schema_version = .ARTIFACT_QUERY_CURSOR_VERSION,
            query_id = specification$query_id,
            artifact_id = source$ref$artifact_id,
            sort_id = sort$sort_id,
            direction = sort$direction,
            values = lapply(statements$cursor_aliases, function(alias) {
                unname(result[[alias]][[nrow(result)]])
            })
        )
    }
    data <- artifactQueryDecodeResult(
        result, projections, columns, statements$status_aliases
    )
    result_bytes <- as.numeric(utils::object.size(data))
    if (nrow(data) > bounds$limit || result_bytes > bounds$max_bytes) {
        artifactQueryAbort(
            "artifact page result exceeds its materialization bounds",
            "multischolar_artifact_query_result_limit_exceeded",
            result_bytes = result_bytes
        )
    }
    artifactQueryAssertRss(bounds$policy, "after_page_query")
    structure(
        list(
            data = data,
            next_cursor = next_cursor,
            has_more = has_more,
            limit = bounds$limit
        ),
        class = c("MultiScholaRArtifactQueryPage", "list")
    )
}
