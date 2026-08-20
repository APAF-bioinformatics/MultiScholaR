makeProjectRegistryFixture <- function(project_root = NULL, policy = NULL) {
    if (is.null(project_root)) {
        project_root <- tempfile("project-registry-fixture-")
        dir.create(project_root, recursive = TRUE)
    }
    identity <- list(
        omic_type = "proteomics",
        omic_label = "proteomics",
        workflow_slug = "diann"
    )
    paths <- buildArtifactPaths(project_root, identity)
    list(
        root = project_root,
        paths = paths,
        registry = newProjectRegistry(
            paths,
            project_id = "project-001",
            resource_policy = policy
        )
    )
}

skipProjectRegistryDependencies <- function() {
    for (package in PROJECT_REGISTRY_REQUIRED_PACKAGES) {
        testthat::skip_if_not_installed(package)
    }
}

projectRegistryTestDigest <- function(letter = "a") {
    paste(rep(letter, 64L), collapse = "")
}

projectRegistryTestTime <- function() {
    artifactRefUtcNow()
}

projectRegistryWriteWorkflow <- function(session, workflow_id = "workflow-001") {
    now <- projectRegistryTestTime()
    projectRegistryWrite(session, "workflow", list(
        workflow_id = workflow_id,
        omic_type = "proteomics",
        omic_label = "proteomics",
        workflow_slug = paste0("diann-", digest::digest(workflow_id, algo = "xxhash32")),
        status = "active",
        created_at = now,
        updated_at = now
    ))
}

projectRegistryWriteRun <- function(
    session,
    workflow_id = "workflow-001",
    run_id = "run-001"
) {
    now <- projectRegistryTestTime()
    projectRegistryWrite(session, "run", list(
        workflow_id = workflow_id,
        run_id = run_id,
        status = "completed",
        action_id = paste0("action-", run_id),
        started_at = now,
        completed_at = now,
        created_at = now,
        updated_at = now
    ))
}

projectRegistryWriteArtifact <- function(
    session,
    artifact_id,
    hydration_ordinal,
    workflow_id = "workflow-001",
    run_id = "run-001",
    generation_id = "generation-001"
) {
    now <- projectRegistryTestTime()
    projectRegistryWrite(session, "artifact", list(
        workflow_id = workflow_id,
        artifact_id = artifact_id,
        run_id = run_id,
        generation_id = generation_id,
        stage_id = "peptide_qc",
        state_role = paste0("role-", hydration_ordinal),
        hydration_ordinal = hydration_ordinal,
        relative_path = file.path(
            "data", "proteomics", "diann", "tables", paste0(artifact_id, ".parquet")
        ),
        codec_id = "data_frame",
        codec_version = 1L,
        payload_schema_id = "test.table",
        payload_schema_version = 1L,
        semantic_digest = projectRegistryTestDigest("a"),
        byte_digest = projectRegistryTestDigest("b"),
        row_count = 10L,
        column_count = 3L,
        payload_bytes = 1024L,
        status = "committed",
        created_at = now,
        updated_at = now
    ))
}

createProjectRegistryMetadataFixture <- function(registry, version, interrupted = FALSE) {
    projectRegistryEnsureParent(registry, "database")
    driver <- duckdb::duckdb(
        dbdir = projectRegistryPath(registry, "database"),
        shared_home = FALSE,
        allow_extensions = FALSE,
        environment_scan = FALSE
    )
    connection <- DBI::dbConnect(driver)
    on.exit(DBI::dbDisconnect(connection, shutdown = TRUE), add = TRUE)
    DBI::dbExecute(connection, paste(
        "CREATE TABLE registry_metadata (",
        "singleton INTEGER PRIMARY KEY, schema_id VARCHAR NOT NULL,",
        "schema_version INTEGER NOT NULL, project_id VARCHAR NOT NULL,",
        "created_at VARCHAR NOT NULL, updated_at VARCHAR NOT NULL)"
    ))
    now <- projectRegistryTestTime()
    projectRegistryExecuteBound(
        connection,
        "INSERT INTO registry_metadata VALUES (?, ?, ?, ?, ?, ?)",
        list(1L, PROJECT_REGISTRY_SCHEMA_ID, version, registry$project_id, now, now)
    )
    if (isTRUE(interrupted)) {
        DBI::dbExecute(connection, paste(
            "CREATE TABLE registry_migrations (",
            "migration_version INTEGER PRIMARY KEY, migration_name VARCHAR NOT NULL,",
            "migration_checksum VARCHAR NOT NULL, applied_at VARCHAR NOT NULL,",
            "package_version VARCHAR NOT NULL)"
        ))
        projectRegistryExecuteBound(
            connection,
            "INSERT INTO registry_migrations VALUES (1, ?, ?, ?, ?)",
            list("interrupted", projectRegistryTestDigest("c"), now, "test")
        )
    }
    invisible(registry)
}

test_that("project registry construction and memory contexts are inert", {
    root <- withr::local_tempdir()
    fixture <- makeProjectRegistryFixture(root)
    expect_s3_class(fixture$registry, "MultiScholaRProjectRegistry")
    expect_length(list.files(root, all.files = TRUE, no.. = TRUE), 0L)
    expect_error(
        projectRegistryAssertLocalFilesystem("nfs4"),
        class = "multischolar_unsupported_registry_filesystem"
    )
    expect_error(
        newProjectRegistry(
            fixture$paths,
            "project-001",
            list(threads = 1000L)
        ),
        class = "multischolar_invalid_registry_resource_policy"
    )

    identity <- list(
        project_id = "project-memory",
        omic_type = "proteomics",
        omic_label = "proteomics",
        workflow_id = "proteomics.gui",
        workflow_type = "DIA",
        workflow_slug = "prot_dia",
        input_format = "diann",
        data_level = "peptide",
        acquisition_mode = "dia"
    )
    context <- WorkflowContext$new(root, identity[c(
        "project_id", "omic_type", "omic_label", "workflow_id"
    )])
    context$bind(identity, list(
        requested_backend = "memory",
        effective_backend = "memory",
        effective_rollout = "memory",
        capability_id = "test.memory",
        capability_version = "1",
        reason_code = "explicit_memory",
        project_state = "new"
    ))
    expect_null(projectRegistryForContext(context))
    expect_length(list.files(root, all.files = TRUE, no.. = TRUE), 0L)

    availability <- stats::setNames(
        c(TRUE, FALSE, TRUE),
        PROJECT_REGISTRY_REQUIRED_PACKAGES
    )
    expect_error(
        projectRegistryAssertDependencies(availability),
        "artifact mode requires optional packages: duckdb",
        class = "multischolar_missing_artifact_dependencies"
    )
})

test_that("fresh registry initialization records bounded and locked resources", {
    skipProjectRegistryDependencies()
    fixture <- makeProjectRegistryFixture()
    session <- initializeProjectRegistry(fixture$registry)
    on.exit(closeProjectRegistry(session), add = TRUE)

    expect_setequal(
        DBI::dbListTables(projectRegistrySessionConnection(session)),
        projectRegistryExpectedTables()
    )
    schema <- projectRegistryQuery(session, "schema")
    expect_identical(schema$schema_id, PROJECT_REGISTRY_SCHEMA_ID)
    expect_identical(as.integer(schema$schema_version), 1L)
    expect_null(session$backup)
    expect_true(file.exists(projectRegistryPath(fixture$registry, "owner")))

    settings <- projectRegistryQuery(
        session,
        "resource_settings",
        list(session_id = session$session_id)
    )
    expect_setequal(
        settings$setting_name,
        c(
            PROJECT_REGISTRY_DUCKDB_SETTINGS,
            "process_rss_limit", "max_result_rows", "max_result_bytes",
            "filesystem_type"
        )
    )
    expect_identical(
        settings$effective_value[settings$setting_name == "temp_directory"],
        fixture$registry$relative_paths$temporary
    )
    expect_false(any(grepl(fixture$root, settings$effective_value, fixed = TRUE)))
    effective <- projectRegistryEffectiveSettings(
        projectRegistrySessionConnection(session)
    )
    expect_silent(projectRegistryValidateEffectiveSettings(effective, fixture$registry))
    expect_error(
        DBI::dbExecute(projectRegistrySessionConnection(session), "SET threads = 64"),
        "configuration has been locked"
    )
    expect_error(
        DBI::dbExecute(projectRegistrySessionConnection(session), "INSTALL httpfs"),
        "extension loading.*disabled"
    )
    inside <- file.path(fixture$root, "inside.csv")
    outside <- tempfile("outside-project-", fileext = ".csv")
    writeLines(c("value", "1"), inside)
    writeLines(c("value", "2"), outside)
    expect_equal(
        projectRegistryFetchBound(
            projectRegistrySessionConnection(session),
            "SELECT * FROM read_csv_auto(?)",
            list(inside)
        )$value,
        1L
    )
    expect_error(
        projectRegistryFetchBound(
            projectRegistrySessionConnection(session),
            "SELECT * FROM read_csv_auto(?)",
            list(outside)
        ),
        "file system operations are disabled"
    )

    expect_true(closeProjectRegistry(session))
    expect_false(file.exists(projectRegistryPath(fixture$registry, "owner")))
    expect_false(closeProjectRegistry(session))
    read_only <- openProjectRegistryReadOnly(fixture$registry)
    on.exit(closeProjectRegistry(read_only), add = TRUE)
    expect_equal(nrow(projectRegistryQuery(read_only, "schema")), 1L)
    expect_error(
        projectRegistryWriteWorkflow(read_only),
        class = "multischolar_registry_read_only"
    )
})

test_that("registry writes enforce complete provenance and relationship integrity", {
    skipProjectRegistryDependencies()
    fixture <- makeProjectRegistryFixture()
    session <- initializeProjectRegistry(fixture$registry)
    on.exit(closeProjectRegistry(session), add = TRUE)
    now <- projectRegistryTestTime()

    projectRegistryWriteWorkflow(session)
    projectRegistryWriteRun(session)
    projectRegistryWrite(session, "source", list(
        workflow_id = "workflow-001",
        run_id = "run-001",
        source_id = "source-001",
        source_role = "primary_input",
        source_uri = "report.tsv'; DROP TABLE projects; --",
        source_digest = projectRegistryTestDigest("c"),
        source_size_bytes = 4096L,
        parser_id = "diann",
        parser_version = "2.0",
        format_id = "tsv",
        data_level = "peptide",
        recorded_at = now
    ))
    projectRegistryWrite(session, "parameter", list(
        workflow_id = "workflow-001",
        run_id = "run-001",
        parameter_id = "parameter-001",
        parameter_name = "q_value",
        value_json = "{\"value\":0.01}",
        value_digest = projectRegistryTestDigest("d"),
        recorded_at = now
    ))
    projectRegistryWrite(session, "software", list(
        workflow_id = "workflow-001",
        run_id = "run-001",
        software_id = "software-001",
        software_name = "MultiScholaR",
        software_version = "0.5.0",
        software_source = "R package",
        software_digest = projectRegistryTestDigest("e"),
        recorded_at = now
    ))
    projectRegistryWriteArtifact(session, "artifact-001", 0L)
    projectRegistryWriteArtifact(session, "artifact-002", 1L)
    projectRegistryWrite(session, "state", list(
        workflow_id = "workflow-001",
        generation_id = "generation-001",
        logical_name = "peptide_qc",
        manifest_relative_path = file.path(
            "state", "proteomics", "diann", "generations", "generation-001.json"
        ),
        status = "current",
        created_at = now,
        updated_at = now
    ))
    projectRegistryWrite(session, "dependency", list(
        workflow_id = "workflow-001",
        artifact_id = "artifact-002",
        depends_on_artifact_id = "artifact-001",
        relationship_type = "derived_from",
        ordinal = 0L,
        recorded_at = now
    ))
    for (direction in c("input", "output")) {
        projectRegistryWrite(session, "run_artifact", list(
            workflow_id = "workflow-001",
            run_id = "run-001",
            artifact_id = if (direction == "input") "artifact-001" else "artifact-002",
            direction = direction,
            artifact_role = "peptide_table",
            ordinal = 0L,
            recorded_at = now
        ))
    }
    projectRegistryWrite(session, "event", list(
        workflow_id = "workflow-001",
        event_id = "event-001",
        generation_id = "generation-001",
        run_id = "run-001",
        event_type = "state_committed",
        event_status = "completed",
        details_json = "{\"rows\":10}",
        recorded_at = now
    ))
    projectRegistryWrite(session, "figure", list(
        workflow_id = "workflow-001",
        figure_id = "figure-001",
        generation_id = "generation-001",
        artifact_id = "artifact-002",
        figure_role = "qc_plot",
        relative_path = "data/proteomics/diann/figures/qc.png",
        content_digest = projectRegistryTestDigest("f"),
        recorded_at = now
    ))
    projectRegistryWrite(session, "metric", list(
        workflow_id = "workflow-001",
        metric_id = "metric-001",
        generation_id = "generation-001",
        run_id = "run-001",
        metric_name = "retained_peptides",
        numeric_value = 10,
        value_json = "{\"count\":10}",
        unit = "peptides",
        recorded_at = now
    ))
    projectRegistryWrite(session, "revision", list(
        workflow_id = "workflow-001",
        revision_id = "revision-001",
        generation_id = "generation-001",
        action_id = "action-run-001",
        revision_status = "committed",
        details_json = "{\"parent\":null}",
        recorded_at = now
    ))

    expect_equal(nrow(projectRegistryQuery(session, "sources", list(
        workflow_id = "workflow-001", run_id = "run-001"
    ))), 1L)
    expect_equal(nrow(projectRegistryQuery(session, "parameters", list(
        workflow_id = "workflow-001", run_id = "run-001"
    ))), 1L)
    expect_equal(nrow(projectRegistryQuery(session, "software", list(
        workflow_id = "workflow-001", run_id = "run-001"
    ))), 1L)
    expect_equal(nrow(projectRegistryQuery(session, "artifacts", list(
        workflow_id = "workflow-001"
    ))), 2L)
    directions <- projectRegistryQuery(session, "run_artifacts", list(
        workflow_id = "workflow-001", run_id = "run-001"
    ))$direction
    expect_setequal(directions, c("input", "output"))
    expect_equal(nrow(projectRegistryQuery(session, "dependencies", list(
        workflow_id = "workflow-001", artifact_id = "artifact-002"
    ))), 1L)
    expect_equal(nrow(projectRegistryQuery(session, "events", list(
        workflow_id = "workflow-001"
    ))), 1L)
    expect_equal(nrow(projectRegistryQuery(session, "figures", list(
        workflow_id = "workflow-001"
    ))), 1L)
    expect_equal(nrow(projectRegistryQuery(session, "metrics", list(
        workflow_id = "workflow-001"
    ))), 1L)
    expect_equal(nrow(projectRegistryQuery(session, "revisions", list(
        workflow_id = "workflow-001"
    ))), 1L)

    expect_error(
        projectRegistryWriteRun(session, workflow_id = "missing-workflow"),
        class = "multischolar_registry_write_failed"
    )
    expect_error(
        projectRegistryWriteArtifact(session, "artifact-003", 1L),
        class = "multischolar_registry_write_failed"
    )
    expect_error(
        projectRegistryWrite(session, "event", list(
            workflow_id = "workflow-001",
            event_id = "event-002",
            event_type = "invalid",
            event_status = "failed",
            details_json = "not-json",
            recorded_at = now
        )),
        class = "multischolar_invalid_registry_json"
    )
})

test_that("registry operations reject SQL, identifier, path, and resource injection", {
    skipProjectRegistryDependencies()
    fixture <- makeProjectRegistryFixture(policy = list(
        max_result_rows = 2L,
        max_result_bytes = 1024^2
    ))
    session <- initializeProjectRegistry(fixture$registry)
    on.exit(closeProjectRegistry(session), add = TRUE)

    expect_error(
        projectRegistryQuery(session, "schema; PRAGMA threads=99"),
        class = "multischolar_unknown_registry_operation"
    )
    expect_error(
        projectRegistryQuery(session, "schema", list(`project_id; DROP` = "x")),
        class = "multischolar_invalid_registry_query"
    )
    expect_error(
        projectRegistryQuery(session, "schema", limit = 1000000000L),
        class = "multischolar_invalid_registry_resource_policy"
    )
    expect_error(
        projectRegistryWrite(session, "artifact", list(
            workflow_id = "workflow-001",
            artifact_id = "artifact-escape",
            generation_id = "generation-escape",
            stage_id = "stage",
            state_role = "role",
            hydration_ordinal = 0L,
            relative_path = "../outside.parquet",
            codec_id = "table",
            codec_version = 1L,
            payload_schema_id = "table",
            payload_schema_version = 1L,
            semantic_digest = projectRegistryTestDigest("a"),
            byte_digest = projectRegistryTestDigest("b"),
            payload_bytes = 1L,
            status = "committed",
            created_at = projectRegistryTestTime(),
            updated_at = projectRegistryTestTime()
        )),
        class = "multischolar_invalid_relative_artifact_path"
    )

    for (index in 1:3) {
        projectRegistryWriteWorkflow(session, paste0("workflow-", index))
    }
    expect_equal(nrow(projectRegistryQuery(session, "workflows")), 2L)
    injected_id <- "workflow-'); DROP TABLE projects; --"
    expect_silent(projectRegistryWriteWorkflow(session, injected_id))
    expect_equal(
        nrow(projectRegistryQuery(
            session,
            "workflows",
            list(workflow_id = injected_id)
        )),
        1L
    )
    expect_true("projects" %in% DBI::dbListTables(
        projectRegistrySessionConnection(session)
    ))
})

test_that("registry migrations back up prior schemas and reject unsafe fixtures", {
    skipProjectRegistryDependencies()

    prior <- makeProjectRegistryFixture()
    createProjectRegistryMetadataFixture(prior$registry, 0L)
    prior_digest <- artifactByteDigest(projectRegistryPath(prior$registry, "database"))
    migrated <- initializeProjectRegistry(prior$registry)
    expect_false(is.null(migrated$backup))
    expect_identical(migrated$backup$sha256, prior_digest)
    expect_true(file.exists(artifactResolveContainedPath(
        prior$root,
        migrated$backup$relative_path,
        must_exist = TRUE
    )))
    expect_equal(projectRegistryQuery(migrated, "schema")$schema_version, 1)
    closeProjectRegistry(migrated)

    future <- makeProjectRegistryFixture()
    createProjectRegistryMetadataFixture(future$registry, 99L)
    future_path <- projectRegistryPath(future$registry, "database", must_exist = TRUE)
    future_digest <- artifactByteDigest(future_path)
    expect_error(
        initializeProjectRegistry(future$registry),
        class = "multischolar_future_registry_schema"
    )
    expect_identical(artifactByteDigest(future_path), future_digest)
    expect_false(dir.exists(projectRegistryPath(future$registry, "backups")))

    malformed <- makeProjectRegistryFixture()
    projectRegistryEnsureParent(malformed$registry, "database")
    malformed_path <- projectRegistryPath(malformed$registry, "database")
    writeBin(charToRaw("not a DuckDB database"), malformed_path)
    malformed_digest <- artifactByteDigest(malformed_path)
    expect_error(
        initializeProjectRegistry(malformed$registry),
        class = "multischolar_registry_open_failed"
    )
    expect_identical(artifactByteDigest(malformed_path), malformed_digest)

    interrupted <- makeProjectRegistryFixture()
    createProjectRegistryMetadataFixture(interrupted$registry, 0L, interrupted = TRUE)
    expect_error(initializeProjectRegistry(interrupted$registry))
    inspected <- projectRegistryInspectExisting(interrupted$registry)
    expect_identical(inspected$version, 0L)
    handle <- projectRegistryConnect(interrupted$registry, TRUE, inspection = TRUE)
    expect_false("projects" %in% DBI::dbListTables(handle$connection))
    projectRegistryDisconnect(handle)
})

test_that("writer ownership is token checked and never recovered by age or PID", {
    skipProjectRegistryDependencies()

    mismatch <- makeProjectRegistryFixture()
    writer <- projectRegistryAcquireWriter(mismatch$registry)
    altered <- writer$owner
    altered$owner_token <- artifactOpaqueId("owner")
    projectRegistryWriteOwner(mismatch$registry, altered)
    expect_error(
        projectRegistryReleaseWriter(writer, writer$owner$owner_token),
        class = "multischolar_registry_owner_token_mismatch"
    )
    projectRegistryWriteOwner(mismatch$registry, writer$owner)
    expect_true(projectRegistryReleaseWriter(writer, writer$owner$owner_token))

    stale <- makeProjectRegistryFixture()
    projectRegistryEnsureParent(stale$registry, "owner")
    stale_owner <- newProjectRegistryOwner(stale$registry)
    stale_owner$pid <- as.integer(Sys.getpid())
    stale_owner$owner_token <- artifactOpaqueId("owner")
    stale_owner$acquired_at <- "2000-01-01T00:00:00.000000Z"
    projectRegistryWriteOwner(stale$registry, stale_owner)
    recovered <- projectRegistryAcquireWriter(stale$registry)
    expect_false(identical(recovered$owner$owner_token, stale_owner$owner_token))
    expect_true(projectRegistryReleaseWriter(
        recovered,
        recovered$owner$owner_token
    ))

    remote <- makeProjectRegistryFixture()
    projectRegistryEnsureParent(remote$registry, "owner")
    remote_owner <- newProjectRegistryOwner(remote$registry)
    remote_owner$host <- paste0(projectRegistryHost(), "-remote")
    remote_owner$acquired_at <- "2000-01-01T00:00:00.000000Z"
    projectRegistryWriteOwner(remote$registry, remote_owner)
    expect_error(
        projectRegistryAcquireWriter(remote$registry),
        class = "multischolar_cross_host_registry_owner"
    )
    unlink(projectRegistryPath(remote$registry, "owner"), force = TRUE)
})

test_that("separate processes exclude live writers and recover abrupt exits", {
    skipProjectRegistryDependencies()
    skip_if_not_installed("processx")
    fixture <- makeProjectRegistryFixture()
    repository <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
    script <- tempfile("project-registry-child-", fileext = ".R")
    writeLines(c(
        "args <- commandArgs(trailingOnly = TRUE)",
        "devtools::load_all(args[[1L]], quiet = TRUE)",
        "identity <- list(",
        "    omic_type = 'proteomics',",
        "    omic_label = 'proteomics',",
        "    workflow_slug = 'diann'",
        ")",
        "paths <- buildArtifactPaths(args[[2L]], identity)",
        "registry <- newProjectRegistry(paths, 'project-001')",
        "session <- initializeProjectRegistry(registry)",
        "writeLines('ready', args[[3L]])",
        "if (identical(args[[4L]], 'abrupt')) quit(save = 'no', status = 0, runLast = FALSE)",
        "while (!file.exists(args[[5L]])) Sys.sleep(0.05)",
        "closeProjectRegistry(session)"
    ), script)

    run_child <- function(mode) {
        ready <- tempfile("project-registry-ready-")
        release <- tempfile("project-registry-release-")
        process <- processx::process$new(
            file.path(R.home("bin"), "Rscript"),
            c("--vanilla", script, repository, fixture$root, ready, mode, release),
            stdout = "|",
            stderr = "2>&1"
        )
        for (attempt in seq_len(600L)) {
            if (file.exists(ready) || !process$is_alive()) break
            Sys.sleep(0.05)
        }
        if (!file.exists(ready)) {
            output <- paste(process$read_all_output_lines(), collapse = "\n")
            testthat::fail(paste("registry child did not become ready:", output))
        }
        list(process = process, ready = ready, release = release)
    }

    live <- run_child("normal")
    expect_error(
        initializeProjectRegistry(fixture$registry),
        class = "multischolar_registry_writer_busy"
    )
    writeLines("release", live$release)
    live$process$wait(timeout = 30000L)
    expect_false(live$process$is_alive())
    expect_identical(live$process$get_exit_status(), 0L)
    after_release <- initializeProjectRegistry(fixture$registry)
    expect_true(closeProjectRegistry(after_release))

    abrupt <- run_child("abrupt")
    abrupt$process$wait(timeout = 30000L)
    expect_false(abrupt$process$is_alive())
    expect_identical(abrupt$process$get_exit_status(), 0L)
    expect_true(file.exists(projectRegistryPath(fixture$registry, "owner")))
    after_abrupt <- initializeProjectRegistry(fixture$registry)
    expect_true(closeProjectRegistry(after_abrupt))
    expect_false(file.exists(projectRegistryPath(fixture$registry, "owner")))
})

test_that("migration interruption rolls back every newly created table", {
    skipProjectRegistryDependencies()
    fixture <- makeProjectRegistryFixture()
    fail_after_schema <- function(stage, migration) {
        expect_identical(stage, "after_schema")
        expect_identical(migration$version, 1L)
        rlang::abort("injected migration interruption", class = "registry_test_interrupt")
    }
    expect_error(
        initializeProjectRegistry(
            fixture$registry,
            failure_injector = fail_after_schema
        ),
        class = "registry_test_interrupt"
    )
    handle <- projectRegistryConnect(fixture$registry, TRUE, inspection = TRUE)
    expect_length(DBI::dbListTables(handle$connection), 0L)
    projectRegistryDisconnect(handle)
    recovered <- initializeProjectRegistry(fixture$registry)
    expect_equal(projectRegistryQuery(recovered, "schema")$schema_version, 1)
    expect_true(closeProjectRegistry(recovered))
})

test_that("migration checksums and independent materialization bounds are enforced", {
    skipProjectRegistryDependencies()

    checksum <- makeProjectRegistryFixture()
    initialized <- initializeProjectRegistry(checksum$registry)
    closeProjectRegistry(initialized)
    driver <- duckdb::duckdb(
        dbdir = projectRegistryPath(checksum$registry, "database", must_exist = TRUE),
        shared_home = FALSE,
        allow_extensions = FALSE,
        environment_scan = FALSE
    )
    connection <- DBI::dbConnect(driver)
    DBI::dbExecute(
        connection,
        "UPDATE registry_migrations SET migration_checksum = 'tampered'"
    )
    DBI::dbDisconnect(connection, shutdown = TRUE)
    database <- projectRegistryPath(checksum$registry, "database", must_exist = TRUE)
    tampered_digest <- artifactByteDigest(database)
    expect_error(
        initializeProjectRegistry(checksum$registry),
        class = "multischolar_malformed_registry_migration_history"
    )
    expect_identical(artifactByteDigest(database), tampered_digest)
    expect_false(dir.exists(projectRegistryPath(checksum$registry, "backups")))

    bounded <- makeProjectRegistryFixture(policy = list(max_result_bytes = 256L))
    bounded_session <- initializeProjectRegistry(bounded$registry)
    expect_error(
        projectRegistryQuery(bounded_session, "resource_settings"),
        class = "multischolar_registry_result_limit_exceeded"
    )
    closeProjectRegistry(bounded_session)

    rss <- makeProjectRegistryFixture(policy = list(
        duckdb_memory_limit_bytes = 1024^2,
        process_rss_limit_bytes = 2 * 1024^2,
        max_result_bytes = 512 * 1024
    ))
    expect_error(
        projectRegistryAssertRss(rss$registry, "test"),
        class = "multischolar_registry_rss_limit_exceeded"
    )
    expect_false(file.exists(projectRegistryPath(rss$registry, "database")))
})
