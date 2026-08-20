# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

projectRegistryExecuteBound <- function(connection, statement, values = list()) {
    result <- DBI::dbSendStatement(connection, statement)
    on.exit(DBI::dbClearResult(result), add = TRUE)
    if (length(values) > 0L) DBI::dbBind(result, values)
    DBI::dbGetRowsAffected(result)
}

projectRegistryFetchBound <- function(connection, statement, values = list(), n = -1L) {
    result <- DBI::dbSendQuery(connection, statement)
    on.exit(DBI::dbClearResult(result), add = TRUE)
    if (length(values) > 0L) DBI::dbBind(result, values)
    DBI::dbFetch(result, n = n)
}

projectRegistryMigrationCatalog <- function() {
    statements <- projectRegistrySchemaStatements()
    statements <- sub(
        "^CREATE TABLE ",
        "CREATE TABLE IF NOT EXISTS ",
        statements
    )
    statements <- sub(
        "^CREATE INDEX ",
        "CREATE INDEX IF NOT EXISTS ",
        statements
    )
    list(list(
        version = 1L,
        name = "initial_project_registry",
        statements = statements,
        checksum = digest::digest(
            paste(statements, collapse = "\n-- statement --\n"),
            algo = "sha256",
            serialize = FALSE
        )
    ))
}

projectRegistryExpectedTables <- function() {
    sort(names(projectRegistryExpectedColumns()))
}

projectRegistryExpectedColumns <- function() {
    specifications <- projectRegistryWriteSpecifications()
    columns <- list(
        registry_metadata = c(
            "singleton", "schema_id", "schema_version", "project_id",
            "created_at", "updated_at"
        ),
        registry_migrations = c(
            "migration_version", "migration_name", "migration_checksum",
            "applied_at", "package_version"
        ),
        registry_resource_settings = c(
            "session_id", "setting_name", "requested_value", "effective_value",
            "unit", "recorded_at"
        ),
        projects = c(
            "project_id", "root_locator", "status", "created_at", "updated_at"
        )
    )
    for (specification in specifications) {
        columns[[specification$table]] <- specification$columns
    }
    columns
}

projectRegistryValidateTableColumns <- function(connection) {
    expected <- projectRegistryExpectedColumns()
    malformed <- names(expected)[vapply(names(expected), function(table) {
        !identical(DBI::dbListFields(connection, table), expected[[table]])
    }, logical(1))]
    if (length(malformed) > 0L) {
        projectRegistryAbort(
            "project registry table structure does not match its migration",
            "multischolar_malformed_registry_schema",
            malformed_tables = malformed
        )
    }
    invisible(TRUE)
}

projectRegistryReadMetadata <- function(connection) {
    projectRegistryFetchBound(
        connection,
        paste0(
            "SELECT singleton, schema_id, schema_version, project_id, ",
            "created_at, updated_at FROM registry_metadata ORDER BY singleton"
        )
    )
}

projectRegistryInspectSchema <- function(connection, project_id) {
    tables <- sort(DBI::dbListTables(connection))
    has_metadata <- "registry_metadata" %in% tables
    if (!has_metadata) {
        if (length(tables) > 0L) {
            projectRegistryAbort(
                "project registry database has tables but no schema metadata",
                "multischolar_malformed_registry_schema",
                tables = tables
            )
        }
        return(list(
            version = 0L,
            created_at = NULL,
            is_fresh = TRUE,
            tables = tables
        ))
    }
    metadata <- tryCatch(
        projectRegistryReadMetadata(connection),
        error = \(error) projectRegistryAbort(
            "project registry schema metadata cannot be read",
            "multischolar_malformed_registry_schema",
            parent = error
        )
    )
    valid <- nrow(metadata) == 1L && identical(as.integer(metadata$singleton), 1L) &&
        identical(metadata$schema_id[[1L]], PROJECT_REGISTRY_SCHEMA_ID) &&
        length(metadata$schema_version) == 1L &&
        !is.na(metadata$schema_version[[1L]]) &&
        metadata$schema_version[[1L]] >= 0L &&
        metadata$schema_version[[1L]] == floor(metadata$schema_version[[1L]]) &&
        identical(metadata$project_id[[1L]], project_id) &&
        artifactRefValidUtc(metadata$created_at[[1L]]) &&
        artifactRefValidUtc(metadata$updated_at[[1L]])
    if (!isTRUE(valid)) {
        projectRegistryAbort(
            "project registry schema metadata is malformed or belongs to another project",
            "multischolar_malformed_registry_schema"
        )
    }
    version <- as.integer(metadata$schema_version[[1L]])
    if (version > PROJECT_REGISTRY_SCHEMA_VERSION) {
        projectRegistryAbort(
            sprintf("project registry schema version %d is newer than this package", version),
            "multischolar_future_registry_schema",
            schema_version = version
        )
    }
    if (version > 0L) {
        projectRegistryValidateMigrationHistory(connection, version)
        missing <- setdiff(projectRegistryExpectedTables(), tables)
        if (length(missing) > 0L) {
            projectRegistryAbort(
                "project registry schema is incomplete",
                "multischolar_malformed_registry_schema",
                missing_tables = missing
            )
        }
        projectRegistryValidateTableColumns(connection)
    }
    list(
        version = version,
        created_at = metadata$created_at[[1L]],
        is_fresh = FALSE,
        tables = tables
    )
}

projectRegistryValidateMigrationHistory <- function(connection, schema_version) {
    if (!"registry_migrations" %in% DBI::dbListTables(connection)) {
        projectRegistryAbort(
            "project registry migration history is missing",
            "multischolar_malformed_registry_migration_history"
        )
    }
    applied <- tryCatch(
        projectRegistryFetchBound(
            connection,
            paste0(
                "SELECT migration_version, migration_name, migration_checksum ",
                "FROM registry_migrations ORDER BY migration_version"
            )
        ),
        error = \(error) projectRegistryAbort(
            "project registry migration history cannot be read",
            "multischolar_malformed_registry_migration_history",
            parent = error
        )
    )
    catalog <- projectRegistryMigrationCatalog()
    expected <- catalog[vapply(catalog, \(item) item$version <= schema_version, logical(1))]
    valid <- nrow(applied) == length(expected)
    if (isTRUE(valid)) {
        for (index in seq_along(expected)) {
            item <- expected[[index]]
            valid <- valid && identical(
                as.integer(applied$migration_version[[index]]),
                item$version
            ) && identical(applied$migration_name[[index]], item$name) &&
                identical(applied$migration_checksum[[index]], item$checksum)
        }
    }
    if (!isTRUE(valid)) {
        projectRegistryAbort(
            "project registry migration history checksum or ordering is invalid",
            "multischolar_malformed_registry_migration_history"
        )
    }
    invisible(applied)
}

projectRegistryBackup <- function(registry, schema_version) {
    registry <- validateProjectRegistry(registry)
    database <- projectRegistryPath(registry, "database", must_exist = TRUE)
    projectRegistryEnsureDirectory(registry, "backups")
    source_digest <- artifactByteDigest(database)
    timestamp <- gsub(
        "[.]",
        "",
        format(Sys.time(), "%Y%m%dT%H%M%OS6", tz = "UTC")
    )
    filename <- sprintf(
        "project-registry-v%d-%s-%s.duckdb",
        schema_version,
        timestamp,
        substr(source_digest, 1L, 12L)
    )
    relative_path <- artifactNormalizeRelativePath(
        file.path(registry$relative_paths$backups, filename)
    )
    path <- artifactResolveContainedPath(registry$project_root, relative_path)
    if (file.exists(path) || !file.copy(database, path, overwrite = FALSE)) {
        projectRegistryAbort(
            "project registry backup could not be created before migration",
            "multischolar_registry_backup_failed"
        )
    }
    backup_digest <- artifactByteDigest(path)
    if (!identical(backup_digest, source_digest)) {
        projectRegistryAbort(
            "project registry backup digest does not match its source",
            "multischolar_registry_backup_failed"
        )
    }
    list(
        relative_path = relative_path,
        schema_version = as.integer(schema_version),
        sha256 = backup_digest,
        created_at = artifactRefUtcNow()
    )
}

projectRegistryInsertMetadata <- function(connection, registry, schema_info, migration) {
    now <- artifactRefUtcNow()
    created_at <- schema_info$created_at
    if (is.null(created_at)) created_at <- now
    projectRegistryExecuteBound(connection, "DELETE FROM registry_metadata")
    projectRegistryExecuteBound(
        connection,
        paste0(
            "INSERT INTO registry_metadata ",
            "(singleton, schema_id, schema_version, project_id, created_at, updated_at) ",
            "VALUES (?, ?, ?, ?, ?, ?)"
        ),
        list(
            1L,
            PROJECT_REGISTRY_SCHEMA_ID,
            migration$version,
            registry$project_id,
            created_at,
            now
        )
    )
    projectRegistryExecuteBound(
        connection,
        paste0(
            "INSERT INTO registry_migrations ",
            "(migration_version, migration_name, migration_checksum, applied_at, ",
            "package_version) VALUES (?, ?, ?, ?, ?)"
        ),
        list(
            migration$version,
            migration$name,
            migration$checksum,
            now,
            projectRegistryPackageVersion()
        )
    )
    projectRegistryExecuteBound(
        connection,
        paste0(
            "INSERT INTO projects ",
            "(project_id, root_locator, status, created_at, updated_at) ",
            "VALUES (?, '.', 'active', ?, ?) ON CONFLICT DO NOTHING"
        ),
        list(registry$project_id, created_at, now)
    )
    invisible(now)
}

projectRegistryMigrate <- function(
    connection,
    registry,
    schema_info,
    failure_injector = NULL
) {
    if (schema_info$version == PROJECT_REGISTRY_SCHEMA_VERSION) {
        return(invisible(FALSE))
    }
    if (!identical(schema_info$version, 0L)) {
        projectRegistryAbort(
            "project registry schema has no supported migration path",
            "multischolar_unsupported_registry_schema"
        )
    }
    migration <- projectRegistryMigrationCatalog()[[1L]]
    DBI::dbBegin(connection)
    committed <- FALSE
    on.exit({
        if (!committed) try(DBI::dbRollback(connection), silent = TRUE)
    }, add = TRUE)
    for (statement in migration$statements) {
        projectRegistryExecuteBound(connection, statement)
    }
    if (!is.null(failure_injector)) {
        if (!is.function(failure_injector)) {
            projectRegistryAbort(
                "project registry migration failure injector must be a function",
                "multischolar_invalid_registry_failure_injector"
            )
        }
        failure_injector("after_schema", migration)
    }
    projectRegistryInsertMetadata(connection, registry, schema_info, migration)
    DBI::dbCommit(connection)
    committed <- TRUE
    invisible(TRUE)
}
