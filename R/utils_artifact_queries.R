# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

artifactQueryAbort <- function(message, class, ...) {
    rlang::abort(
        message,
        class = c(class, "multischolar_artifact_query_error"),
        ...
    )
}

artifactQuerySpecification <- function(descriptor, operation_id) {
    descriptor <- validateArtifactWorkflowDescriptor(descriptor)
    if (!workflowCapabilityScalarString(operation_id) ||
        !operation_id %in% names(descriptor$queries)) {
        artifactQueryAbort(
            "artifact query operation is not owned by the workflow descriptor",
            "multischolar_unknown_artifact_query_operation"
        )
    }
    descriptor$queries[[operation_id]]
}

artifactQueryValidateReference <- function(store, ref, state_role, descriptor) {
    store <- validateArtifactStore(store)
    descriptor <- validateArtifactWorkflowDescriptor(descriptor)
    ref <- artifactStoreNormalizeRef(ref)
    owned <- identical(ref$logical_key$state_role, state_role) &&
        identical(ref$logical_key$omic_type, descriptor$identity$omic_type) &&
        identical(
            ref$logical_key$workflow_slug,
            descriptor$identity$workflow_slug
        )
    if (!isTRUE(owned)) {
        artifactQueryAbort(
            "artifact query operation does not own this workflow or state role",
            "multischolar_artifact_query_role_mismatch"
        )
    }
    pin_path <- artifactWorkflowStateDescriptorPinPath(store)
    if (!file.exists(artifactStoreResolveFile(store, pin_path))) {
        artifactQueryAbort(
            "artifact query workflow has no immutable descriptor pin",
            "multischolar_artifact_query_descriptor_pin_missing"
        )
    }
    pin <- artifactStoreReadJson(store, pin_path)
    expected_contract <- list(
        descriptor_id = descriptor$descriptor_id,
        descriptor_version = descriptor$descriptor_version,
        descriptor_digest = descriptor$descriptor_digest
    )
    pin_version <- workflowStateVersionValue(pin$schema_version)
    pin_valid <- identical(pin$schema, ARTIFACT_DESCRIPTOR_PIN_SCHEMA) &&
        identical(pin_version, ARTIFACT_DESCRIPTOR_PIN_VERSION) &&
        identical(pin$project_id, store$project_id) &&
        identical(pin$workflow_id, descriptor$identity$workflow_id) &&
        identical(pin$contract, expected_contract)
    if (!isTRUE(pin_valid)) {
        artifactQueryAbort(
            "artifact query descriptor differs from the persisted workflow pin",
            "multischolar_artifact_query_descriptor_pin_mismatch"
        )
    }
    sidecar_path <- artifactStoreManagedPaths(
        store, ref$logical_key, ref$artifact_id
    )$sidecar
    sidecar <- artifactStoreReadSidecar(store, sidecar_path, validate_payload = FALSE)
    if (!identical(sidecar$artifact_ref, ref)) {
        artifactQueryAbort(
            "artifact query reference does not match its immutable sidecar",
            "multischolar_artifact_query_ref_mismatch"
        )
    }
    payload_path <- artifactStoreResolveFile(
        store, ref$relative_path, must_exist = TRUE
    )
    bytes <- unname(as.numeric(file.info(payload_path)$size))
    if (!identical(bytes, as.numeric(ref$shape$bytes))) {
        artifactQueryAbort(
            "artifact query payload size differs from its immutable reference",
            "multischolar_artifact_query_payload_mismatch"
        )
    }
    list(ref = ref, sidecar = sidecar, payload_path = payload_path)
}

artifactQueryColumnMap <- function(codec_metadata) {
    valid <- is.list(codec_metadata) &&
        identical(codec_metadata$kind, "data.frame") &&
        is.character(codec_metadata$logical_names) &&
        is.list(codec_metadata$columns) &&
        length(codec_metadata$logical_names) == length(codec_metadata$columns) &&
        anyDuplicated(codec_metadata$logical_names) == 0L &&
        .artifactRowOrderColumn %in% codec_metadata$internal_columns
    if (!isTRUE(valid)) {
        artifactQueryAbort(
            "artifact query requires table codec metadata with unique logical columns",
            "multischolar_unsupported_artifact_query_payload"
        )
    }
    physical <- vapply(
        codec_metadata$columns, `[[`, character(1), "physical_name"
    )
    types <- vapply(codec_metadata$columns, `[[`, character(1), "r_type")
    descriptors <- stats::setNames(
        codec_metadata$columns, codec_metadata$logical_names
    )
    list(
        physical = stats::setNames(physical, codec_metadata$logical_names),
        types = stats::setNames(types, codec_metadata$logical_names),
        descriptors = descriptors
    )
}

artifactQueryNormalizeProjections <- function(specification, projections, columns) {
    if (is.null(projections)) projections <- specification$projections
    valid <- is.character(projections) && length(projections) > 0L &&
        !anyNA(projections) && all(nzchar(projections)) &&
        anyDuplicated(projections) == 0L &&
        all(projections %in% specification$projections) &&
        all(projections %in% names(columns$physical))
    if (!isTRUE(valid)) {
        artifactQueryAbort(
            "artifact query projection is absent from the descriptor or payload",
            "multischolar_invalid_artifact_query_projection"
        )
    }
    unname(projections)
}

artifactQueryValidateValue <- function(value, type, vector = FALSE) {
    length_valid <- if (isTRUE(vector)) length(value) >= 1L else length(value) == 1L
    type_valid <- switch(type,
        character = is.character(value),
        integer = is.integer(value),
        double = is.numeric(value),
        logical = is.logical(value),
        FALSE
    )
    valid <- length_valid && type_valid && !anyNA(value)
    if (identical(type, "double")) valid <- valid && all(is.finite(value))
    if (!isTRUE(valid)) {
        artifactQueryAbort(
            "artifact query filter value does not match its declared type",
            "multischolar_invalid_artifact_query_filter"
        )
    }
    unname(value)
}

artifactQueryFilterSql <- function(
    specification,
    filters,
    columns,
    connection
) {
    if (!is.list(filters) ||
        (length(filters) > 0L && is.null(names(filters))) ||
        any(!nzchar(names(filters))) || anyDuplicated(names(filters)) > 0L ||
        !all(names(filters) %in% names(specification$filters))) {
        artifactQueryAbort(
            "artifact query filters are not descriptor-owned named operations",
            "multischolar_invalid_artifact_query_filter"
        )
    }
    clauses <- character()
    values <- list()
    quote <- function(value) {
        as.character(DBI::dbQuoteIdentifier(connection, value))
    }
    for (filter_id in names(filters)) {
        declaration <- specification$filters[[filter_id]]
        request <- filters[[filter_id]]
        if (!is.list(request) ||
            !identical(names(request), c("operator", "value")) ||
            !workflowCapabilityScalarString(request$operator) ||
            !request$operator %in% declaration$operators) {
            artifactQueryAbort(
                sprintf("artifact query filter '%s' is malformed", filter_id),
                "multischolar_invalid_artifact_query_filter"
            )
        }
        column <- columns$physical[[declaration$column]]
        actual_type <- columns$types[[declaration$column]]
        compatible_type <- identical(actual_type, declaration$type) ||
            (actual_type %in% c("factor", "ordered") &&
                identical(declaration$type, "character"))
        if (is.null(column) || !isTRUE(compatible_type)) {
            artifactQueryAbort(
                sprintf("artifact query filter '%s' disagrees with payload type", filter_id),
                "multischolar_artifact_query_schema_mismatch"
            )
        }
        operator <- request$operator
        if (identical(operator, "is_missing")) {
            value <- artifactQueryValidateValue(request$value, "logical")
            status_column <- columns$descriptors[[declaration$column]]$status_name
            clause <- if (is.null(status_column)) {
                paste0(
                    quote(column),
                    if (isTRUE(value)) " IS NULL" else " IS NOT NULL"
                )
            } else {
                paste0(
                    quote(status_column),
                    if (isTRUE(value)) " = 1" else " <> 1"
                )
            }
            clauses <- c(clauses, clause)
            next
        }
        vector <- operator %in% c("in", "between")
        value <- artifactQueryValidateValue(
            request$value, declaration$type, vector = vector
        )
        if (identical(operator, "between") && length(value) != 2L) {
            artifactQueryAbort(
                "artifact between filter requires exactly two values",
                "multischolar_invalid_artifact_query_filter"
            )
        }
        if (identical(operator, "in") && length(value) > 1000L) {
            artifactQueryAbort(
                "artifact membership filter exceeds its value bound",
                "multischolar_artifact_query_filter_limit_exceeded"
            )
        }
        sql <- switch(operator,
            equal = paste0(quote(column), " = ?"),
            gte = paste0(quote(column), " >= ?"),
            lte = paste0(quote(column), " <= ?"),
            between = paste0(quote(column), " BETWEEN ? AND ?"),
            `in` = paste0(
                quote(column), " IN (", paste(rep("?", length(value)), collapse = ", "),
                ")"
            )
        )
        clauses <- c(clauses, sql)
        values <- c(values, as.list(value))
    }
    list(
        sql = if (length(clauses) == 0L) {
            ""
        } else {
            paste0(" WHERE ", paste(clauses, collapse = " AND "))
        },
        values = values
    )
}

artifactQueryNormalizeBounds <- function(specification, limit, resource_policy) {
    policy <- normalizeProjectRegistryPolicy(resource_policy)
    maximum_rows <- min(specification$max_rows, policy$max_result_rows)
    maximum_bytes <- min(specification$max_bytes, policy$max_result_bytes)
    if (is.null(limit)) limit <- maximum_rows
    valid <- length(limit) == 1L && is.numeric(limit) && !is.na(limit) &&
        is.finite(limit) && limit >= 1L && limit == floor(limit) &&
        limit <= maximum_rows
    if (!isTRUE(valid)) {
        artifactQueryAbort(
            "artifact query row limit exceeds its descriptor or resource bound",
            "multischolar_artifact_query_row_limit_exceeded"
        )
    }
    list(
        limit = as.integer(limit),
        max_bytes = as.numeric(maximum_bytes),
        policy = policy
    )
}

artifactQueryAssertRss <- function(policy, stage) {
    rss <- projectRegistryProcessRss()
    if (!is.na(rss) && rss > policy$process_rss_limit_bytes) {
        artifactQueryAbort(
            sprintf("process RSS exceeds the artifact query limit at '%s'", stage),
            "multischolar_artifact_query_rss_limit_exceeded",
            rss_bytes = rss,
            limit_bytes = policy$process_rss_limit_bytes
        )
    }
    invisible(rss)
}

artifactQueryConnect <- function(store, policy) {
    projectRegistryRequireDependencies()
    temporary <- artifactNormalizeRelativePath(store$relative_paths$duckdb_tmp)
    artifactStoreEnsureDirectory(store, temporary)
    config <- list(
        threads = as.character(policy$threads),
        memory_limit = projectRegistryByteSetting(policy$duckdb_memory_limit_bytes),
        max_temp_directory_size = projectRegistryByteSetting(
            policy$temp_size_limit_bytes
        ),
        temp_directory = artifactStoreResolveFile(store, temporary),
        preserve_insertion_order = "false",
        enable_external_file_cache = "false",
        autoinstall_known_extensions = "false",
        autoload_known_extensions = "false",
        allow_community_extensions = "false",
        allow_unsigned_extensions = "false"
    )
    driver <- tryCatch(
        duckdb::duckdb(
            dbdir = ":memory:",
            config = config,
            shared_home = FALSE,
            allow_extensions = FALSE,
            environment_scan = FALSE
        ),
        error = function(error) artifactQueryAbort(
            "bounded artifact query could not initialize DuckDB",
            "multischolar_artifact_query_open_failed",
            parent = error
        )
    )
    connection <- tryCatch(
        DBI::dbConnect(driver),
        error = function(error) {
            try(duckdb::duckdb_shutdown(driver), silent = TRUE)
            artifactQueryAbort(
                "bounded artifact query could not connect to DuckDB",
                "multischolar_artifact_query_open_failed",
                parent = error
            )
        }
    )
    tryCatch(
        {
            DBI::dbExecute(
                connection,
                "SET allowed_directories = ?",
                params = list(list(store$project_root))
            )
            DBI::dbExecute(connection, "SET lock_configuration = true")
        },
        error = function(error) {
            try(DBI::dbDisconnect(connection, shutdown = TRUE), silent = TRUE)
            artifactQueryAbort(
                "bounded artifact query security policy could not be locked",
                "multischolar_artifact_query_security_failed",
                parent = error
            )
        }
    )
    list(driver = driver, connection = connection)
}

artifactQueryDisconnect <- function(handle) {
    if (!is.null(handle$connection) && DBI::dbIsValid(handle$connection)) {
        DBI::dbDisconnect(handle$connection, shutdown = TRUE)
    } else if (!is.null(handle$driver)) {
        try(duckdb::duckdb_shutdown(handle$driver), silent = TRUE)
    }
    invisible(TRUE)
}

artifactQueryStatements <- function(
    connection,
    projections,
    order_by,
    columns,
    filter_sql
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
        sprintf(".multischolar_query_status_%04d", seq_along(status_projections)),
        status_projections
    )
    status_terms <- vapply(status_projections, function(projection) {
        paste0(
            quote(columns$descriptors[[projection]]$status_name),
            " AS ", quote(status_aliases[[projection]])
        )
    }, character(1))
    select <- paste(
        c(projection_terms, status_terms),
        collapse = ", "
    )
    order <- paste(
        quote(c(columns$physical[order_by], .artifactRowOrderColumn)),
        collapse = ", "
    )
    inner <- paste0(
        "SELECT ", select, " FROM read_parquet(?)", filter_sql,
        " ORDER BY ", order, " LIMIT ?"
    )
    byte_terms <- paste0(
        "octet_length(encode(COALESCE(CAST(",
        quote(c(projections, unname(status_aliases))),
        " AS VARCHAR), '')))"
    )
    estimate <- paste0(
        "WITH bounded_artifact_query AS (", inner, ") ",
        "SELECT COUNT(*) AS row_count, COALESCE(SUM(",
        paste(byte_terms, collapse = " + "),
        "), 0) AS payload_bytes FROM bounded_artifact_query"
    )
    list(
        query = inner,
        estimate = estimate,
        status_aliases = status_aliases
    )
}

artifactQueryDecodeResult <- function(result, projections, columns, status_aliases) {
    decoded <- lapply(projections, function(projection) {
        descriptor <- columns$descriptors[[projection]]
        payload <- list()
        payload[[descriptor$physical_name]] <- result[[projection]]
        if (!is.null(descriptor$status_name)) {
            payload[[descriptor$status_name]] <- result[[status_aliases[[projection]]]]
        }
        artifactDecodeColumn(payload, descriptor)
    })
    names(decoded) <- projections
    as.data.frame(decoded, optional = TRUE, stringsAsFactors = FALSE)
}

queryArtifactRef <- function(
    store,
    ref,
    descriptor,
    operation_id,
    projections = NULL,
    filters = list(),
    limit = NULL,
    resource_policy = NULL
) {
    specification <- artifactQuerySpecification(descriptor, operation_id)
    source <- artifactQueryValidateReference(
        store, ref, specification$state_role, descriptor
    )
    columns <- artifactQueryColumnMap(source$sidecar$codec_metadata)
    projections <- artifactQueryNormalizeProjections(
        specification, projections, columns
    )
    if (!all(specification$order_by %in% names(columns$physical))) {
        artifactQueryAbort(
            "artifact query ordering columns are absent from the payload",
            "multischolar_artifact_query_schema_mismatch"
        )
    }
    bounds <- artifactQueryNormalizeBounds(specification, limit, resource_policy)
    artifactQueryAssertRss(bounds$policy, "before_query")
    handle <- artifactQueryConnect(store, bounds$policy)
    on.exit(artifactQueryDisconnect(handle), add = TRUE)
    filter <- artifactQueryFilterSql(
        specification, filters, columns, handle$connection
    )
    statements <- artifactQueryStatements(
        handle$connection,
        projections,
        specification$order_by,
        columns,
        filter$sql
    )
    values <- c(list(source$payload_path), filter$values, list(bounds$limit))
    estimate <- tryCatch(
        DBI::dbGetQuery(handle$connection, statements$estimate, params = values),
        error = function(error) artifactQueryAbort(
            "bounded artifact query byte estimation failed",
            "multischolar_artifact_query_failed",
            parent = error
        )
    )
    estimated_bytes <- as.numeric(estimate$payload_bytes[[1L]]) * 2 +
        as.numeric(estimate$row_count[[1L]]) * 256
    if (estimated_bytes > bounds$max_bytes) {
        artifactQueryAbort(
            "artifact query exceeds its pre-materialization byte limit",
            "multischolar_artifact_query_byte_limit_exceeded",
            estimated_bytes = estimated_bytes
        )
    }
    result <- tryCatch(
        DBI::dbGetQuery(handle$connection, statements$query, params = values),
        error = function(error) artifactQueryAbort(
            "bounded artifact query failed",
            "multischolar_artifact_query_failed",
            parent = error
        )
    )
    result <- artifactQueryDecodeResult(
        result, projections, columns, statements$status_aliases
    )
    result_bytes <- as.numeric(utils::object.size(result))
    if (nrow(result) > bounds$limit || result_bytes > bounds$max_bytes) {
        artifactQueryAbort(
            "artifact query result exceeds its materialization bounds",
            "multischolar_artifact_query_result_limit_exceeded",
            result_bytes = result_bytes
        )
    }
    artifactQueryAssertRss(bounds$policy, "after_query")
    result
}
