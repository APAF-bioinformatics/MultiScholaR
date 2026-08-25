# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

projectRegistryWriteSpecifications <- function() {
    list(
        workflow = list(
            table = "workflows",
            columns = c(
                "project_id", "workflow_id", "omic_type", "omic_label",
                "workflow_slug", "status", "created_at", "updated_at"
            ),
            optional = character()
        ),
        run = list(
            table = "workflow_runs",
            columns = c(
                "project_id", "workflow_id", "run_id", "status", "action_id",
                "started_at", "completed_at", "created_at", "updated_at"
            ),
            optional = c("action_id", "completed_at")
        ),
        source = list(
            table = "provenance_sources",
            columns = c(
                "project_id", "workflow_id", "run_id", "source_id",
                "source_role", "source_uri", "source_digest", "source_size_bytes",
                "parser_id", "parser_version", "format_id", "data_level",
                "recorded_at"
            ),
            optional = c("source_uri", "source_size_bytes")
        ),
        parameter = list(
            table = "provenance_parameters",
            columns = c(
                "project_id", "workflow_id", "run_id", "parameter_id",
                "parameter_name", "value_json", "value_digest", "recorded_at"
            ),
            optional = character()
        ),
        software = list(
            table = "provenance_software",
            columns = c(
                "project_id", "workflow_id", "run_id", "software_id",
                "software_name", "software_version", "software_source",
                "software_digest", "recorded_at"
            ),
            optional = "software_digest"
        ),
        artifact = list(
            table = "artifacts",
            columns = c(
                "project_id", "workflow_id", "artifact_id", "run_id",
                "generation_id", "stage_id", "state_role", "hydration_ordinal",
                "relative_path", "codec_id", "codec_version", "payload_schema_id",
                "payload_schema_version", "semantic_digest", "byte_digest",
                "row_count", "column_count", "payload_bytes", "status",
                "created_at", "updated_at"
            ),
            optional = c("run_id", "row_count", "column_count")
        ),
        state = list(
            table = "workflow_states",
            columns = c(
                "project_id", "workflow_id", "generation_id",
                "parent_generation_id", "logical_name", "manifest_relative_path",
                "status", "created_at", "updated_at"
            ),
            optional = "parent_generation_id"
        ),
        dependency = list(
            table = "artifact_dependencies",
            columns = c(
                "project_id", "workflow_id", "artifact_id",
                "depends_on_artifact_id", "relationship_type", "ordinal",
                "recorded_at"
            ),
            optional = character()
        ),
        run_artifact = list(
            table = "run_artifacts",
            columns = c(
                "project_id", "workflow_id", "run_id", "artifact_id",
                "direction", "artifact_role", "ordinal", "recorded_at"
            ),
            optional = character()
        ),
        event = list(
            table = "workflow_events",
            columns = c(
                "project_id", "workflow_id", "event_id", "generation_id",
                "run_id", "event_type", "event_status", "details_json",
                "recorded_at"
            ),
            optional = c("generation_id", "run_id")
        ),
        figure = list(
            table = "workflow_figures",
            columns = c(
                "project_id", "workflow_id", "figure_id", "generation_id",
                "artifact_id", "figure_role", "relative_path", "content_digest",
                "recorded_at"
            ),
            optional = "artifact_id"
        ),
        metric = list(
            table = "workflow_metrics",
            columns = c(
                "project_id", "workflow_id", "metric_id", "generation_id",
                "run_id", "metric_name", "numeric_value", "value_json", "unit",
                "recorded_at"
            ),
            optional = c("generation_id", "run_id", "numeric_value", "value_json", "unit")
        ),
        revision = list(
            table = "registry_revisions",
            columns = c(
                "project_id", "workflow_id", "revision_id", "generation_id",
                "action_id", "expected_parent_generation_id", "revision_status",
                "details_json", "recorded_at"
            ),
            optional = c("generation_id", "expected_parent_generation_id")
        )
    )
}

projectRegistryStatusValues <- function(column, operation) {
    key <- paste(operation, column, sep = ":")
    switch(
        key,
        "workflow:status" = c("created", "active", "completed", "failed", "stale", "archived"),
        "run:status" = c("created", "running", "completed", "failed", "cancelled"),
        "artifact:status" = c("staged", "validated", "committed", "trashed"),
        "state:status" = c("staged", "current", "historical", "stale", "failed"),
        "run_artifact:direction" = c("input", "output"),
        NULL
    )
}

projectRegistryNullableValue <- function(column) {
    if (grepl("(_bytes|_count|_ordinal|numeric_value)$", column) ||
        identical(column, "ordinal")) {
        return(NA_real_)
    }
    NA_character_
}

projectRegistryValidateRecordValue <- function(
    registry,
    operation,
    column,
    value,
    nullable
) {
    if (is.null(value) || (length(value) == 1L && is.na(value))) {
        if (isTRUE(nullable)) return(projectRegistryNullableValue(column))
        projectRegistryAbort(
            sprintf("project registry field '%s' is required", column),
            "multischolar_invalid_registry_record"
        )
    }
    if (length(value) != 1L || is.list(value) || is.factor(value)) {
        projectRegistryAbort(
            sprintf("project registry field '%s' must be scalar", column),
            "multischolar_invalid_registry_record"
        )
    }
    if (is.character(value) &&
        nchar(value, type = "bytes") > registry$resource_policy$max_result_bytes) {
        projectRegistryAbort(
            sprintf("project registry field '%s' exceeds its byte limit", column),
            "multischolar_registry_result_limit_exceeded"
        )
    }
    if (grepl("(_at|started_at|completed_at)$", column)) {
        if (!artifactRefValidUtc(value)) {
            projectRegistryAbort(
                sprintf("project registry field '%s' must be a UTC timestamp", column),
                "multischolar_invalid_registry_timestamp"
            )
        }
        return(value)
    }
    if (grepl("(_digest|byte_digest|semantic_digest)$", column)) {
        artifactRefValidateDigest(value, column)
        return(value)
    }
    if (grepl("(_json|value_json|details_json)$", column)) {
        if (!workflowCapabilityScalarString(value) || !jsonlite::validate(value)) {
            projectRegistryAbort(
                sprintf("project registry field '%s' must contain valid JSON", column),
                "multischolar_invalid_registry_json"
            )
        }
        return(value)
    }
    if (grepl("relative_path$", column)) {
        relative <- artifactNormalizeRelativePath(value)
        artifactResolveContainedPath(registry$project_root, relative)
        return(relative)
    }
    if (column %in% c("codec_version", "payload_schema_version")) {
        return(projectRegistryPositiveInteger(value, column))
    }
    if (grepl("(_bytes|_count|_ordinal)$", column) || identical(column, "ordinal")) {
        valid <- is.numeric(value) && is.finite(value) && value >= 0 && value == floor(value)
        if (!isTRUE(valid)) {
            projectRegistryAbort(
                sprintf("project registry field '%s' must be a non-negative integer", column),
                "multischolar_invalid_registry_record"
            )
        }
        return(as.numeric(value))
    }
    if (identical(column, "numeric_value")) {
        if (!is.numeric(value) || !is.finite(value)) {
            projectRegistryAbort(
                "project registry metric value must be finite",
                "multischolar_invalid_registry_record"
            )
        }
        return(as.numeric(value))
    }
    allowed <- projectRegistryStatusValues(column, operation)
    if (!is.null(allowed) && !value %in% allowed) {
        projectRegistryAbort(
            sprintf("project registry field '%s' has an unsupported value", column),
            "multischolar_invalid_registry_status"
        )
    }
    if (!workflowCapabilityScalarString(value)) {
        projectRegistryAbort(
            sprintf("project registry field '%s' must be a non-empty string", column),
            "multischolar_invalid_registry_record"
        )
    }
    value
}

projectRegistryNormalizeRecord <- function(registry, operation, record) {
    specifications <- projectRegistryWriteSpecifications()
    if (!workflowCapabilityScalarString(operation) ||
        !operation %in% names(specifications)) {
        projectRegistryAbort(
            "project registry write operation is not package-owned",
            "multischolar_unknown_registry_operation"
        )
    }
    if (!is.list(record) || is.null(names(record)) ||
        any(!nzchar(names(record))) || anyDuplicated(names(record)) > 0L) {
        projectRegistryAbort(
            "project registry record must be an unambiguous named list",
            "multischolar_invalid_registry_record"
        )
    }
    specification <- specifications[[operation]]
    allowed_input <- setdiff(specification$columns, "project_id")
    if (length(setdiff(names(record), allowed_input)) > 0L) {
        projectRegistryAbort(
            "project registry record contains unsupported fields",
            "multischolar_invalid_registry_record"
        )
    }
    required <- setdiff(allowed_input, specification$optional)
    if (length(setdiff(required, names(record))) > 0L) {
        projectRegistryAbort(
            "project registry record is missing required fields",
            "multischolar_invalid_registry_record"
        )
    }
    values <- c(list(project_id = registry$project_id), record)
    values <- lapply(specification$columns, function(column) {
        value <- values[[column]]
        projectRegistryValidateRecordValue(
            registry,
            operation,
            column,
            value,
            nullable = column %in% specification$optional
        )
    })
    names(values) <- specification$columns
    list(specification = specification, values = values)
}

projectRegistryAssertRelation <- function(connection, table, filters, relationship) {
    quote <- \(values) as.character(DBI::dbQuoteIdentifier(connection, values))
    statement <- paste0(
        "SELECT COUNT(*) AS matches FROM ", quote(table), " WHERE ",
        paste(paste0(quote(names(filters)), " = ?"), collapse = " AND ")
    )
    matches <- projectRegistryFetchBound(
        connection,
        statement,
        unname(filters)
    )$matches[[1L]]
    if (!identical(as.numeric(matches), 1)) {
        projectRegistryAbort(
            sprintf("project registry relationship '%s' does not exist", relationship),
            "multischolar_invalid_registry_relationship",
            relationship = relationship
        )
    }
    invisible(TRUE)
}

projectRegistryValidateRelationships <- function(connection, operation, values) {
    scoped <- \(extra) c(
        list(project_id = values$project_id, workflow_id = values$workflow_id),
        extra
    )
    relation <- function(table, extra, label) {
        projectRegistryAssertRelation(connection, table, scoped(extra), label)
    }
    if (operation == "artifact" && !is.na(values$run_id)) {
        relation("workflow_runs", list(run_id = values$run_id), "artifact_run")
    }
    if (operation == "state" && !is.na(values$parent_generation_id)) {
        relation(
            "workflow_states",
            list(generation_id = values$parent_generation_id),
            "state_parent"
        )
    }
    if (operation %in% c("event", "metric") && !is.na(values$run_id)) {
        relation("workflow_runs", list(run_id = values$run_id), paste0(operation, "_run"))
    }
    if (operation %in% c("event", "figure", "metric", "revision") &&
        !is.na(values$generation_id)) {
        relation(
            "workflow_states",
            list(generation_id = values$generation_id),
            paste0(operation, "_generation")
        )
    }
    if (operation == "figure" && !is.na(values$artifact_id)) {
        relation("artifacts", list(artifact_id = values$artifact_id), "figure_artifact")
    }
    if (operation == "revision" && !is.na(values$expected_parent_generation_id)) {
        relation(
            "workflow_states",
            list(generation_id = values$expected_parent_generation_id),
            "revision_expected_parent"
        )
    }
    invisible(TRUE)
}

projectRegistryWrite <- function(session, operation, record) {
    session <- validateProjectRegistrySession(session, write = TRUE)
    registry <- session$registry
    projectRegistryAssertRss(registry, paste0("before_write_", operation))
    normalized <- projectRegistryNormalizeRecord(registry, operation, record)
    connection <- projectRegistrySessionConnection(session, write = TRUE)
    projectRegistryValidateRelationships(
        connection,
        operation,
        normalized$values
    )
    columns <- as.character(DBI::dbQuoteIdentifier(
        connection,
        normalized$specification$columns
    ))
    table <- as.character(DBI::dbQuoteIdentifier(
        connection,
        normalized$specification$table
    ))
    statement <- paste0(
        "INSERT INTO ", table, " (", paste(columns, collapse = ", "), ") VALUES (",
        paste(rep("?", length(columns)), collapse = ", "), ")"
    )
    affected <- tryCatch(
        projectRegistryExecuteBound(
            connection,
            statement,
            unname(normalized$values)
        ),
        error = \(error) projectRegistryAbort(
            sprintf("project registry '%s' write failed", operation),
            "multischolar_registry_write_failed",
            parent = error
        )
    )
    projectRegistryAssertRss(registry, paste0("after_write_", operation))
    invisible(affected)
}

projectRegistryQuerySpecifications <- function() {
    list(
        schema = list(
            table = "registry_metadata",
            columns = c("schema_id", "schema_version", "project_id", "created_at", "updated_at"),
            filters = "project_id",
            required = character(),
            order = "schema_version"
        ),
        resource_settings = list(
            table = "registry_resource_settings",
            columns = c(
                "session_id", "setting_name", "requested_value", "effective_value",
                "unit", "recorded_at"
            ),
            filters = "session_id",
            required = character(),
            order = c("recorded_at", "setting_name"),
            project_scoped = FALSE
        ),
        workflows = list(
            table = "workflows",
            columns = c(
                "project_id", "workflow_id", "omic_type", "omic_label",
                "workflow_slug", "status", "created_at", "updated_at"
            ),
            filters = c("project_id", "workflow_id", "omic_type", "status"),
            required = character(),
            order = c("omic_type", "workflow_id")
        ),
        runs = list(
            table = "workflow_runs",
            columns = c(
                "project_id", "workflow_id", "run_id", "status", "action_id",
                "started_at", "completed_at", "created_at", "updated_at"
            ),
            filters = c("project_id", "workflow_id", "run_id", "status"),
            required = "workflow_id",
            order = c("started_at", "run_id")
        ),
        sources = list(
            table = "provenance_sources",
            columns = c(
                "project_id", "workflow_id", "run_id", "source_id", "source_role",
                "source_uri", "source_digest", "source_size_bytes", "parser_id",
                "parser_version", "format_id", "data_level", "recorded_at"
            ),
            filters = c("project_id", "workflow_id", "run_id", "source_role"),
            required = c("workflow_id", "run_id"),
            order = c("source_role", "source_id")
        ),
        parameters = list(
            table = "provenance_parameters",
            columns = c(
                "project_id", "workflow_id", "run_id", "parameter_id",
                "parameter_name", "value_json", "value_digest", "recorded_at"
            ),
            filters = c("project_id", "workflow_id", "run_id", "parameter_name"),
            required = c("workflow_id", "run_id"),
            order = c("parameter_name", "parameter_id")
        ),
        software = list(
            table = "provenance_software",
            columns = c(
                "project_id", "workflow_id", "run_id", "software_id",
                "software_name", "software_version", "software_source",
                "software_digest", "recorded_at"
            ),
            filters = c("project_id", "workflow_id", "run_id", "software_name"),
            required = c("workflow_id", "run_id"),
            order = c("software_name", "software_id")
        ),
        artifacts = list(
            table = "artifacts",
            columns = c(
                "project_id", "workflow_id", "artifact_id", "run_id",
                "generation_id", "stage_id", "state_role", "hydration_ordinal",
                "relative_path", "codec_id", "payload_schema_id", "semantic_digest",
                "byte_digest", "payload_bytes", "status", "created_at", "updated_at"
            ),
            filters = c(
                "project_id", "workflow_id", "artifact_id", "generation_id",
                "stage_id", "state_role", "status"
            ),
            required = "workflow_id",
            order = c("generation_id", "hydration_ordinal", "artifact_id")
        ),
        states = list(
            table = "workflow_states",
            columns = c(
                "project_id", "workflow_id", "generation_id", "parent_generation_id",
                "logical_name", "manifest_relative_path", "status", "created_at",
                "updated_at"
            ),
            filters = c("project_id", "workflow_id", "generation_id", "status"),
            required = "workflow_id",
            order = c("created_at", "generation_id")
        ),
        dependencies = list(
            table = "artifact_dependencies",
            columns = c(
                "project_id", "workflow_id", "artifact_id",
                "depends_on_artifact_id", "relationship_type", "ordinal",
                "recorded_at"
            ),
            filters = c("project_id", "workflow_id", "artifact_id"),
            required = c("workflow_id", "artifact_id"),
            order = c("ordinal", "depends_on_artifact_id")
        ),
        run_artifacts = list(
            table = "run_artifacts",
            columns = c(
                "project_id", "workflow_id", "run_id", "artifact_id", "direction",
                "artifact_role", "ordinal", "recorded_at"
            ),
            filters = c("project_id", "workflow_id", "run_id", "direction"),
            required = c("workflow_id", "run_id"),
            order = c("direction", "artifact_role", "ordinal")
        ),
        events = list(
            table = "workflow_events",
            columns = c(
                "project_id", "workflow_id", "event_id", "generation_id", "run_id",
                "event_type", "event_status", "details_json", "recorded_at"
            ),
            filters = c("project_id", "workflow_id", "generation_id", "run_id", "event_type"),
            required = "workflow_id",
            order = c("recorded_at", "event_id")
        ),
        figures = list(
            table = "workflow_figures",
            columns = c(
                "project_id", "workflow_id", "figure_id", "generation_id",
                "artifact_id", "figure_role", "relative_path", "content_digest",
                "recorded_at"
            ),
            filters = c("project_id", "workflow_id", "generation_id", "figure_role"),
            required = "workflow_id",
            order = c("generation_id", "figure_role", "figure_id")
        ),
        metrics = list(
            table = "workflow_metrics",
            columns = c(
                "project_id", "workflow_id", "metric_id", "generation_id", "run_id",
                "metric_name", "numeric_value", "value_json", "unit", "recorded_at"
            ),
            filters = c("project_id", "workflow_id", "generation_id", "run_id", "metric_name"),
            required = "workflow_id",
            order = c("recorded_at", "metric_name", "metric_id")
        ),
        revisions = list(
            table = "registry_revisions",
            columns = c(
                "project_id", "workflow_id", "revision_id", "generation_id",
                "action_id", "expected_parent_generation_id", "revision_status",
                "details_json", "recorded_at"
            ),
            filters = c("project_id", "workflow_id", "generation_id", "action_id"),
            required = "workflow_id",
            order = c("recorded_at", "revision_id")
        )
    )
}

#' Resolve a legacy public workflow filter when it is unambiguous
#'
#' @param session An open project registry session.
#' @param operation Package-owned registry query operation.
#' @param filters Candidate query filters.
#'
#' @return Query filters with an exact registry workflow ID when resolvable.
#' @noRd
projectRegistryResolveWorkflowFilter <- function(session, operation, filters) {
    specification <- projectRegistryQuerySpecifications()[[operation]]
    supports_workflow <- !is.null(specification) &&
        !identical(operation, "workflows") &&
        "workflow_id" %in% specification$filters
    if (!isTRUE(supports_workflow)) return(filters)
    if (!is.list(filters) ||
        !workflowCapabilityScalarString(filters$workflow_id)) {
        return(filters)
    }
    workflow_id <- filters$workflow_id
    connection <- projectRegistrySessionConnection(session)
    exact <- projectRegistryFetchBound(
        connection,
        paste0(
            "SELECT workflow_id FROM workflows ",
            "WHERE project_id = ? AND workflow_id = ? LIMIT 1"
        ),
        list(session$registry$project_id, workflow_id)
    )
    if (nrow(exact) == 1L) return(filters)
    scoped <- projectRegistryFetchBound(
        connection,
        paste0(
            "SELECT workflow_id FROM workflows ",
            "WHERE project_id = ? AND starts_with(workflow_id, ?) ",
            "ORDER BY workflow_id LIMIT 2"
        ),
        list(session$registry$project_id, paste0(workflow_id, "::"))
    )
    if (nrow(scoped) > 1L) {
        projectRegistryAbort(
            "public workflow ID matches multiple registry workflows",
            "multischolar_ambiguous_registry_workflow_id"
        )
    }
    if (nrow(scoped) == 1L) filters$workflow_id <- scoped$workflow_id[[1L]]
    filters
}

projectRegistryNormalizeQuery <- function(session, operation, filters, limit) {
    specifications <- projectRegistryQuerySpecifications()
    if (!workflowCapabilityScalarString(operation) ||
        !operation %in% names(specifications)) {
        projectRegistryAbort(
            "project registry query operation is not package-owned",
            "multischolar_unknown_registry_operation"
        )
    }
    if (!is.list(filters) ||
        (length(filters) > 0L && is.null(names(filters))) ||
        any(!nzchar(names(filters))) || anyDuplicated(names(filters)) > 0L) {
        projectRegistryAbort(
            "project registry query filters must be an unambiguous named list",
            "multischolar_invalid_registry_query"
        )
    }
    specification <- specifications[[operation]]
    project_scoped <- specification$project_scoped
    if (is.null(project_scoped)) project_scoped <- TRUE
    caller_allowed <- setdiff(specification$filters, "project_id")
    if (length(setdiff(names(filters), caller_allowed)) > 0L ||
        length(setdiff(specification$required, names(filters))) > 0L) {
        projectRegistryAbort(
            "project registry query has unsupported or missing filters",
            "multischolar_invalid_registry_query"
        )
    }
    if (isTRUE(project_scoped)) {
        filters <- c(list(project_id = session$registry$project_id), filters)
    }
    valid_filter <- vapply(
        filters,
        workflowCapabilityScalarString,
        logical(1)
    )
    if (!all(valid_filter)) {
        projectRegistryAbort(
            "project registry query filters must be scalar strings",
            "multischolar_invalid_registry_query"
        )
    }
    maximum <- session$registry$resource_policy$max_result_rows
    if (is.null(limit)) limit <- maximum
    limit <- as.integer(projectRegistryPositiveInteger(
        limit,
        "query_limit",
        maximum = maximum
    ))
    list(specification = specification, filters = filters, limit = limit)
}

projectRegistryEstimateQueryBytes <- function(
    connection,
    specification,
    where,
    order,
    values
) {
    columns <- as.character(DBI::dbQuoteIdentifier(
        connection,
        specification$columns
    ))
    table <- as.character(DBI::dbQuoteIdentifier(
        connection,
        specification$table
    ))
    byte_terms <- paste0(
        "octet_length(encode(COALESCE(CAST(",
        columns,
        " AS VARCHAR), '')))"
    )
    statement <- paste0(
        "SELECT COUNT(*) AS row_count, COALESCE(SUM(row_bytes), 0) AS payload_bytes ",
        "FROM (SELECT ", paste(byte_terms, collapse = " + "),
        " AS row_bytes FROM ", table, where, " ORDER BY ", order,
        " LIMIT ?) AS bounded_registry_query"
    )
    estimate <- projectRegistryFetchBound(connection, statement, values)
    as.numeric(estimate$payload_bytes[[1L]]) * 2 +
        as.numeric(estimate$row_count[[1L]]) * 256
}

projectRegistryQuery <- function(session, operation, filters = list(), limit = NULL) {
    session <- validateProjectRegistrySession(session)
    registry <- session$registry
    projectRegistryAssertRss(registry, paste0("before_query_", operation))
    filters <- projectRegistryResolveWorkflowFilter(session, operation, filters)
    normalized <- projectRegistryNormalizeQuery(session, operation, filters, limit)
    connection <- projectRegistrySessionConnection(session)
    specification <- normalized$specification
    quote <- \(values) as.character(DBI::dbQuoteIdentifier(connection, values))
    where <- if (length(normalized$filters) == 0L) {
        ""
    } else {
        paste0(
            " WHERE ",
            paste(paste0(quote(names(normalized$filters)), " = ?"), collapse = " AND ")
        )
    }
    statement <- paste0(
        "SELECT ", paste(quote(specification$columns), collapse = ", "),
        " FROM ", quote(specification$table), where,
        " ORDER BY ", paste(quote(specification$order), collapse = ", "),
        " LIMIT ?"
    )
    values <- c(unname(normalized$filters), list(normalized$limit))
    order <- paste(quote(specification$order), collapse = ", ")
    estimated_bytes <- projectRegistryEstimateQueryBytes(
        connection,
        specification,
        where,
        order,
        values
    )
    if (estimated_bytes > registry$resource_policy$max_result_bytes) {
        projectRegistryAbort(
            "project registry query exceeds its pre-materialization byte limit",
            "multischolar_registry_result_limit_exceeded",
            estimated_bytes = estimated_bytes
        )
    }
    result <- tryCatch(
        projectRegistryFetchBound(connection, statement, values),
        error = \(error) projectRegistryAbort(
            sprintf("project registry '%s' query failed", operation),
            "multischolar_registry_query_failed",
            parent = error
        )
    )
    result_bytes <- as.numeric(utils::object.size(result))
    if (result_bytes > registry$resource_policy$max_result_bytes) {
        projectRegistryAbort(
            "project registry query result exceeds its materialization limit",
            "multischolar_registry_result_limit_exceeded",
            result_bytes = result_bytes
        )
    }
    projectRegistryAssertRss(registry, paste0("after_query_", operation))
    result
}
