# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

PROJECT_REGISTRY_DUCKDB_SETTINGS <- c(
    "allow_community_extensions",
    "allow_unsigned_extensions",
    "allowed_directories",
    "autoload_known_extensions",
    "autoinstall_known_extensions",
    "enable_external_access",
    "enable_external_file_cache",
    "lock_configuration",
    "max_temp_directory_size",
    "memory_limit",
    "preserve_insertion_order",
    "temp_directory",
    "threads"
)

projectRegistryByteSetting <- function(value) {
    paste0(format(value, scientific = FALSE, trim = TRUE), "B")
}

projectRegistryConnectionConfig <- function(registry, inspection = FALSE) {
    registry <- validateProjectRegistry(registry)
    policy <- registry$resource_policy
    config <- list(
        threads = as.character(policy$threads),
        memory_limit = projectRegistryByteSetting(policy$duckdb_memory_limit_bytes),
        max_temp_directory_size = projectRegistryByteSetting(
            policy$temp_size_limit_bytes
        ),
        preserve_insertion_order = if (policy$preserve_insertion_order) {
            "true"
        } else {
            "false"
        },
        enable_external_file_cache = "false"
    )
    if (!isTRUE(inspection)) {
        config$temp_directory <- projectRegistryPath(registry, "temporary")
    }
    config$autoinstall_known_extensions <- "false"
    config$autoload_known_extensions <- "false"
    config$allow_community_extensions <- "false"
    config$allow_unsigned_extensions <- "false"
    config
}

projectRegistryConfigureSecurity <- function(connection, registry) {
    projectRegistryExecuteBound(
        connection,
        "SET allowed_directories = ?",
        list(list(registry$project_root))
    )
    projectRegistryExecuteBound(
        connection,
        "SET enable_external_access = false"
    )
    projectRegistryExecuteBound(
        connection,
        "SET lock_configuration = true"
    )
    invisible(TRUE)
}

projectRegistryConnect <- function(registry, read_only, inspection = FALSE) {
    registry <- validateProjectRegistry(registry)
    projectRegistryRequireDependencies()
    database <- projectRegistryPath(
        registry,
        "database",
        must_exist = isTRUE(read_only)
    )
    driver <- tryCatch(
        duckdb::duckdb(
            dbdir = database,
            read_only = isTRUE(read_only),
            config = projectRegistryConnectionConfig(registry, inspection),
            shared_home = FALSE,
            allow_extensions = FALSE,
            environment_scan = FALSE
        ),
        error = \(error) projectRegistryAbort(
            "project registry DuckDB driver could not be opened",
            "multischolar_registry_open_failed",
            parent = error
        )
    )
    connection <- tryCatch(
        DBI::dbConnect(driver),
        error = \(error) {
            try(duckdb::duckdb_shutdown(driver), silent = TRUE)
            projectRegistryAbort(
                "project registry DuckDB connection could not be opened",
                "multischolar_registry_open_failed",
                parent = error
            )
        }
    )
    tryCatch(
        projectRegistryConfigureSecurity(connection, registry),
        error = \(error) {
            try(DBI::dbDisconnect(connection, shutdown = TRUE), silent = TRUE)
            projectRegistryAbort(
                "project registry DuckDB security configuration failed",
                "multischolar_registry_security_configuration_failed",
                parent = error
            )
        }
    )
    list(driver = driver, connection = connection)
}

projectRegistryDisconnect <- function(handle) {
    if (!is.list(handle) ||
        !all(c("driver", "connection") %in% names(handle))) {
        projectRegistryAbort(
            "project registry connection handle is invalid",
            "multischolar_invalid_registry_connection"
        )
    }
    if (DBI::dbIsValid(handle$connection)) {
        DBI::dbDisconnect(handle$connection, shutdown = TRUE)
    } else {
        try(duckdb::duckdb_shutdown(handle$driver), silent = TRUE)
    }
    invisible(TRUE)
}

projectRegistryEffectiveSettings <- function(connection) {
    placeholders <- paste(rep("?", length(PROJECT_REGISTRY_DUCKDB_SETTINGS)), collapse = ", ")
    projectRegistryFetchBound(
        connection,
        paste0(
            "SELECT name, value FROM duckdb_settings() WHERE name IN (",
            placeholders,
            ") ORDER BY name"
        ),
        as.list(PROJECT_REGISTRY_DUCKDB_SETTINGS)
    )
}

projectRegistryParseByteSetting <- function(value) {
    matched <- regexec(
        "^([0-9]+(?:[.][0-9]+)?) (B|KiB|MiB|GiB|TiB)$",
        value
    )
    parts <- regmatches(value, matched)[[1L]]
    if (length(parts) != 3L) return(NA_real_)
    multiplier <- switch(
        parts[[3L]],
        B = 1,
        KiB = 1024,
        MiB = 1024^2,
        GiB = 1024^3,
        TiB = 1024^4
    )
    as.numeric(parts[[2L]]) * multiplier
}

projectRegistryValidateEffectiveSettings <- function(settings, registry = NULL) {
    missing <- setdiff(PROJECT_REGISTRY_DUCKDB_SETTINGS, settings$name)
    if (length(missing) > 0L || anyDuplicated(settings$name) > 0L) {
        projectRegistryAbort(
            "DuckDB did not expose every required project registry setting",
            "multischolar_registry_resource_configuration_failed",
            missing_settings = missing
        )
    }
    values <- stats::setNames(as.character(settings$value), settings$name)
    disabled <- c(
        "allow_community_extensions", "allow_unsigned_extensions",
        "autoload_known_extensions", "autoinstall_known_extensions",
        "enable_external_access", "enable_external_file_cache"
    )
    if (any(tolower(values[disabled]) != "false") ||
        tolower(values[["lock_configuration"]]) != "true") {
        projectRegistryAbort(
            "DuckDB project registry security settings are not locked down",
            "multischolar_registry_security_configuration_failed"
        )
    }
    if (!is.null(registry)) {
        registry <- validateProjectRegistry(registry)
        policy <- registry$resource_policy
        memory_bytes <- projectRegistryParseByteSetting(values[["memory_limit"]])
        temp_bytes <- projectRegistryParseByteSetting(
            values[["max_temp_directory_size"]]
        )
        valid_resources <- !is.na(memory_bytes) &&
            memory_bytes <= policy$duckdb_memory_limit_bytes &&
            !is.na(temp_bytes) && temp_bytes <= policy$temp_size_limit_bytes &&
            as.integer(values[["threads"]]) <= policy$threads &&
            identical(
                tolower(values[["preserve_insertion_order"]]),
                tolower(as.character(policy$preserve_insertion_order))
            ) && identical(
                normalizePath(values[["temp_directory"]], winslash = "/", mustWork = TRUE),
                projectRegistryPath(registry, "temporary", must_exist = TRUE)
            ) && grepl(registry$project_root, values[["allowed_directories"]], fixed = TRUE)
        if (!isTRUE(valid_resources)) {
            projectRegistryAbort(
                "DuckDB project registry resource settings exceed their policy",
                "multischolar_registry_resource_configuration_failed"
            )
        }
    }
    invisible(values)
}

projectRegistryRequestedSetting <- function(registry, setting_name) {
    policy <- registry$resource_policy
    switch(
        setting_name,
        threads = as.character(policy$threads),
        memory_limit = as.character(policy$duckdb_memory_limit_bytes),
        max_temp_directory_size = as.character(policy$temp_size_limit_bytes),
        temp_directory = registry$relative_paths$temporary,
        allowed_directories = ".",
        preserve_insertion_order = as.character(policy$preserve_insertion_order),
        lock_configuration = "true",
        "false"
    )
}

projectRegistrySettingUnit <- function(setting_name) {
    switch(
        setting_name,
        memory_limit = "bytes",
        max_temp_directory_size = "bytes",
        temp_directory = "project_relative_path",
        allowed_directories = "project_root_allowlist",
        threads = "count",
        preserve_insertion_order = "logical",
        lock_configuration = "logical",
        "logical"
    )
}

projectRegistryRecordSettings <- function(connection, registry, session_id) {
    settings <- projectRegistryEffectiveSettings(connection)
    projectRegistryValidateEffectiveSettings(settings, registry)
    now <- artifactRefUtcNow()
    rows <- lapply(seq_len(nrow(settings)), function(index) {
        name <- settings$name[[index]]
        list(
            session_id,
            name,
            projectRegistryRequestedSetting(registry, name),
            if (name %in% c("temp_directory", "allowed_directories")) {
                projectRegistryRequestedSetting(registry, name)
            } else {
                as.character(settings$value[[index]])
            },
            projectRegistrySettingUnit(name),
            now
        )
    })
    rows <- c(rows, list(
        list(
            session_id, "process_rss_limit", 
            as.character(registry$resource_policy$process_rss_limit_bytes),
            as.character(registry$resource_policy$process_rss_limit_bytes),
            "bytes", now
        ),
        list(
            session_id, "max_result_rows",
            as.character(registry$resource_policy$max_result_rows),
            as.character(registry$resource_policy$max_result_rows),
            "rows", now
        ),
        list(
            session_id, "max_result_bytes",
            as.character(registry$resource_policy$max_result_bytes),
            as.character(registry$resource_policy$max_result_bytes),
            "bytes", now
        ),
        list(
            session_id, "filesystem_type", registry$filesystem_type,
            registry$filesystem_type, "filesystem", now
        )
    ))
    DBI::dbBegin(connection)
    committed <- FALSE
    on.exit({
        if (!committed) try(DBI::dbRollback(connection), silent = TRUE)
    }, add = TRUE)
    statement <- paste0(
        "INSERT INTO registry_resource_settings ",
        "(session_id, setting_name, requested_value, effective_value, unit, recorded_at) ",
        "VALUES (?, ?, ?, ?, ?, ?)"
    )
    for (row in rows) projectRegistryExecuteBound(connection, statement, row)
    DBI::dbCommit(connection)
    committed <- TRUE
    invisible(settings)
}

newProjectRegistrySession <- function(
    registry,
    handle,
    read_only,
    session_id,
    writer = NULL,
    backup = NULL
) {
    session <- new.env(parent = emptyenv())
    session$registry <- registry
    session$handle <- handle
    session$read_only <- isTRUE(read_only)
    session$session_id <- session_id
    session$writer <- writer
    session$backup <- backup
    session$closed <- FALSE
    class(session) <- c("MultiScholaRProjectRegistrySession", "environment")
    session
}

validateProjectRegistrySession <- function(session, write = FALSE) {
    valid <- inherits(session, "MultiScholaRProjectRegistrySession") &&
        is.environment(session) && !isTRUE(session$closed) &&
        is.list(session$handle) && DBI::dbIsValid(session$handle$connection)
    if (!isTRUE(valid)) {
        projectRegistryAbort(
            "project registry session is closed or invalid",
            "multischolar_invalid_registry_session"
        )
    }
    if (isTRUE(write) && (isTRUE(session$read_only) || is.null(session$writer))) {
        projectRegistryAbort(
            "project registry write requires the exclusive writer guard",
            "multischolar_registry_read_only"
        )
    }
    session
}

projectRegistrySessionConnection <- function(session, write = FALSE) {
    session <- validateProjectRegistrySession(session, write = write)
    session$handle$connection
}

closeProjectRegistry <- function(session) {
    if (!inherits(session, "MultiScholaRProjectRegistrySession") ||
        !is.environment(session)) {
        projectRegistryAbort(
            "project registry session is invalid",
            "multischolar_invalid_registry_session"
        )
    }
    if (isTRUE(session$closed)) return(invisible(FALSE))
    disconnect_error <- NULL
    if (!is.null(session$handle)) {
        tryCatch(
            projectRegistryDisconnect(session$handle),
            error = \(error) disconnect_error <<- error
        )
        if (is.null(disconnect_error)) session$handle <- NULL
    }
    release_error <- NULL
    if (is.null(disconnect_error) && !is.null(session$writer)) {
        tryCatch(
            projectRegistryReleaseWriter(
                session$writer,
                session$writer$owner$owner_token
            ),
            error = \(error) release_error <<- error
        )
        if (is.null(release_error)) session$writer <- NULL
    }
    session$closed <- is.null(session$handle) && is.null(session$writer)
    if (!is.null(disconnect_error)) stop(disconnect_error)
    if (!is.null(release_error)) stop(release_error)
    invisible(TRUE)
}

projectRegistryInspectExisting <- function(registry) {
    handle <- projectRegistryConnect(registry, read_only = TRUE, inspection = TRUE)
    on.exit(projectRegistryDisconnect(handle), add = TRUE)
    projectRegistryInspectSchema(handle$connection, registry$project_id)
}

initializeProjectRegistry <- function(registry, failure_injector = NULL) {
    registry <- validateProjectRegistry(registry)
    projectRegistryRequireDependencies()
    projectRegistryAssertRss(registry, "before_initialize")
    writer <- projectRegistryAcquireWriter(registry)
    writer_owned <- TRUE
    on.exit({
        if (writer_owned) {
            try(
                projectRegistryReleaseWriter(writer, writer$owner$owner_token),
                silent = TRUE
            )
        }
    }, add = TRUE)
    database_exists <- file.exists(projectRegistryPath(registry, "database"))
    schema_info <- if (database_exists) {
        projectRegistryInspectExisting(registry)
    } else {
        list(version = 0L, created_at = NULL, is_fresh = TRUE, tables = character())
    }
    backup <- NULL
    if (database_exists && schema_info$version < PROJECT_REGISTRY_SCHEMA_VERSION) {
        backup <- projectRegistryBackup(registry, schema_info$version)
    }
    projectRegistryEnsureDirectory(registry, "temporary")
    handle <- projectRegistryConnect(registry, read_only = FALSE)
    handle_owned <- TRUE
    on.exit({
        if (handle_owned) try(projectRegistryDisconnect(handle), silent = TRUE)
    }, add = TRUE)
    projectRegistryMigrate(
        handle$connection,
        registry,
        schema_info,
        failure_injector = failure_injector
    )
    validated <- projectRegistryInspectSchema(handle$connection, registry$project_id)
    if (!identical(validated$version, PROJECT_REGISTRY_SCHEMA_VERSION)) {
        projectRegistryAbort(
            "project registry migration did not reach the current schema",
            "multischolar_registry_migration_failed"
        )
    }
    projectRegistryRecordSettings(
        handle$connection,
        registry,
        writer$owner$session_id
    )
    projectRegistryAssertRss(registry, "after_initialize")
    session <- newProjectRegistrySession(
        registry,
        handle,
        read_only = FALSE,
        session_id = writer$owner$session_id,
        writer = writer,
        backup = backup
    )
    handle_owned <- FALSE
    writer_owned <- FALSE
    session
}

openProjectRegistryReadOnly <- function(registry) {
    registry <- validateProjectRegistry(registry)
    projectRegistryRequireDependencies()
    projectRegistryAssertRss(registry, "before_read_only_open")
    handle <- projectRegistryConnect(registry, read_only = TRUE)
    handle_owned <- TRUE
    on.exit({
        if (handle_owned) try(projectRegistryDisconnect(handle), silent = TRUE)
    }, add = TRUE)
    schema_info <- projectRegistryInspectSchema(handle$connection, registry$project_id)
    if (!identical(schema_info$version, PROJECT_REGISTRY_SCHEMA_VERSION)) {
        projectRegistryAbort(
            "project registry requires migration before read-only use",
            "multischolar_registry_requires_migration"
        )
    }
    settings <- projectRegistryEffectiveSettings(handle$connection)
    projectRegistryValidateEffectiveSettings(settings, registry)
    session <- newProjectRegistrySession(
        registry,
        handle,
        read_only = TRUE,
        session_id = artifactOpaqueId("session")
    )
    handle_owned <- FALSE
    session
}
