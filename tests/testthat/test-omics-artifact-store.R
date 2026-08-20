makeArtifactStoreFixture <- function(project_root, generation_id = "generation-001") {
    identity <- list(
        omic_type = "proteomics",
        omic_label = "proteomics",
        workflow_slug = "diann"
    )
    paths <- buildArtifactPaths(project_root, identity)
    store <- newArtifactStore(paths, project_id = "project-001")
    logical_key <- list(
        project_id = "project-001",
        omic_type = "proteomics",
        workflow_slug = "diann",
        stage_id = "peptide_qc",
        state_role = "filtered_peptides",
        generation_id = generation_id
    )
    value <- data.frame(
        feature_id = c("F1", "F2", "F3"),
        sample = factor(c("S1", "S2", "S1"), levels = c("S2", "S1")),
        value = c(1, NA_real_, Inf),
        stringsAsFactors = FALSE,
        row.names = c("row-c", "row-a", "row-b")
    )
    list(
        paths = paths,
        store = store,
        logical_key = logical_key,
        value = value,
        encoded = encodeArtifactTable(value, stable_key = "feature_id", owner = "fixture")
    )
}

failArtifactStoreAt <- function(target_stage) {
    force(target_stage)
    function(stage, context) {
        if (identical(stage, target_stage)) {
            rlang::abort(
                paste("injected artifact failure at", stage),
                class = "multischolar_test_artifact_failure"
            )
        }
        invisible(context)
    }
}

test_that("ArtifactStore construction is inert and successful writes are exact", {
    project_root <- withr::local_tempdir()
    fixture <- makeArtifactStoreFixture(project_root)
    expect_false(dir.exists(artifactPath(fixture$paths, "tables")))
    expect_false(dir.exists(artifactPath(fixture$paths, "generations")))

    ref <- artifactStoreWriteParquet(
        fixture$store,
        fixture$encoded,
        fixture$logical_key
    )
    managed <- artifactStoreManagedPaths(
        fixture$store,
        fixture$logical_key,
        ref$artifact_id
    )
    payload_path <- artifactStoreResolveFile(
        fixture$store,
        managed$payload,
        must_exist = TRUE
    )
    expect_s3_class(ref, "MultiScholaRArtifactRef")
    expect_identical(ref$status, "committed")
    expect_identical(ref$shape$kind, "data.frame")
    expect_identical(ref$shape$rows, 3L)
    expect_identical(ref$shape$columns, 3L)
    expect_identical(ref$shape$payloads, 1L)
    expect_identical(ref$shape$bytes, unname(as.numeric(file.info(payload_path)$size)))
    expect_identical(ref$hash_policy$byte$digest, artifactByteDigest(payload_path))
    expect_identical(
        ref$hash_policy$semantic$digest,
        fixture$encoded$metadata$semantic_digest
    )

    sidecar <- artifactStoreReadSidecar(
        fixture$store,
        managed$sidecar,
        validate_payload = TRUE
    )
    expect_identical(sidecar$artifact_ref, ref)
    expect_identical(sidecar$registration$logical_key, fixture$logical_key)
    expect_identical(sidecar$codec_metadata, fixture$encoded$metadata)
    expect_match(sidecar$guarantees$visibility, "atomic")
    expect_match(sidecar$guarantees$durability, "not_claimed")
    expect_false(grepl(project_root, paste(capture.output(str(sidecar)), collapse = "")))
    expect_false(file.exists(artifactStoreResolveFile(fixture$store, managed$intent)))

    restored <- decodeArtifactRectangular(
        arrow::read_parquet(payload_path, as_data_frame = FALSE),
        sidecar$codec_metadata
    )
    expect_identical(restored, fixture$value)
    expect_error(
        artifactStoreWriteParquet(
            fixture$store,
            fixture$encoded,
            fixture$logical_key
        ),
        class = "multischolar_duplicate_artifact_generation"
    )
})

test_that("all write boundaries reconcile idempotently", {
    stages <- c(
        "before_write", "after_temp_write", "after_validation",
        "after_payload_rename", "after_sidecar_rename"
    )
    expected_states <- c(
        before_write = "intent_only",
        after_temp_write = "unvalidated_temporary",
        after_validation = "validated_temporary",
        after_payload_rename = "payload_published_sidecar_pending",
        after_sidecar_rename = "committed"
    )
    for (stage in stages) {
        project_root <- withr::local_tempdir()
        fixture <- makeArtifactStoreFixture(project_root)
        expect_error(
            artifactStoreWriteParquet(
                fixture$store,
                fixture$encoded,
                fixture$logical_key,
                failure_injector = failArtifactStoreAt(stage)
            ),
            class = "multischolar_test_artifact_failure",
            info = stage
        )
        initial <- reconcileArtifactStore(fixture$store)
        expect_identical(initial$state, unname(expected_states[[stage]]), info = stage)
        first <- reconcileArtifactStore(fixture$store, repair = TRUE)
        second <- reconcileArtifactStore(fixture$store, repair = TRUE)
        expect_identical(first, second, info = stage)
        if (!identical(stage, "before_write")) {
            expect_identical(first$state, "committed", info = stage)
            expect_identical(first$registry_state, "not_compared", info = stage)
        } else {
            expect_identical(first$state, "intent_only", info = stage)
        }
    }
})

test_that("abrupt child termination leaves a recoverable explicit intent", {
    skip_if(Sys.which("Rscript") == "")
    project_root <- withr::local_tempdir()
    script <- tempfile(fileext = ".R")
    repo_root <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
    writeLines(c(
        "args <- commandArgs(trailingOnly = TRUE)",
        "devtools::load_all(args[[1L]], quiet = TRUE)",
        "identity <- list(",
        "    omic_type = 'proteomics',",
        "    omic_label = 'proteomics',",
        "    workflow_slug = 'diann'",
        ")",
        "paths <- buildArtifactPaths(args[[2L]], identity)",
        "store <- newArtifactStore(paths, 'project-001')",
        "key <- list(",
        "    project_id = 'project-001',",
        "    omic_type = 'proteomics',",
        "    workflow_slug = 'diann',",
        "    stage_id = 'peptide_qc',",
        "    state_role = 'filtered_peptides',",
        "    generation_id = 'generation-abrupt'",
        ")",
        "value <- data.frame(id = c('A', 'B'), value = c(1, 2))",
        "encoded <- encodeArtifactTable(value, stable_key = 'id')",
        "terminate <- function(stage, context) {",
        "    if (stage == 'after_payload_rename') {",
        "        quit(save = 'no', status = 73L, runLast = FALSE)",
        "    }",
        "}",
        "artifactStoreWriteParquet(store, encoded, key, failure_injector = terminate)"
    ), script)
    output <- tempfile(fileext = ".log")
    status <- suppressWarnings(system2(
        Sys.which("Rscript"),
        c("--vanilla", script, repo_root, project_root),
        stdout = output,
        stderr = output
    ))
    expect_identical(as.integer(status), 73L)

    fixture <- makeArtifactStoreFixture(project_root, "generation-abrupt")
    inventory <- reconcileArtifactStore(fixture$store)
    expect_identical(inventory$state, "payload_published_sidecar_pending")
    recovered <- reconcileArtifactStore(fixture$store, repair = TRUE)
    expect_identical(recovered$state, "committed")
})

test_that("write failures retain bounded evidence without permission changes", {
    scenarios <- list(
        disk_full = "No space left on device",
        quota = "Disk quota exceeded",
        permission = "Permission denied"
    )
    for (scenario in names(scenarios)) {
        project_root <- withr::local_tempdir()
        fixture <- makeArtifactStoreFixture(project_root)
        failing_writer <- function(x, sink, ...) {
            if (!identical(scenario, "permission")) {
                connection <- file(sink, open = "wb")
                on.exit(close(connection), add = TRUE)
                writeBin(charToRaw("partial parquet"), connection)
            }
            stop(scenarios[[scenario]], call. = FALSE)
        }
        expect_error(
            artifactStoreWriteParquet(
                fixture$store,
                fixture$encoded,
                fixture$logical_key,
                write_parquet_fn = failing_writer
            ),
            regexp = scenarios[[scenario]],
            class = "multischolar_artifact_write_failed"
        )
        inventory <- reconcileArtifactStore(fixture$store, repair = TRUE)
        expect_true(
            inventory$state %in% c("intent_only", "unvalidated_temporary"),
            info = scenario
        )
        expect_false(dir.exists(file.path(project_root, "unmanaged")))
    }
})

test_that("temporary paths are unique, colocated, and never replace a final", {
    project_root <- withr::local_tempdir()
    first <- makeArtifactStoreFixture(project_root, "generation-temp-1")
    second <- makeArtifactStoreFixture(project_root, "generation-temp-2")
    for (fixture in list(first, second)) {
        expect_error(
            artifactStoreWriteParquet(
                fixture$store,
                fixture$encoded,
                fixture$logical_key,
                failure_injector = failArtifactStoreAt("after_validation")
            ),
            class = "multischolar_test_artifact_failure"
        )
    }
    intent_paths <- artifactStoreIntentPaths(first$store)
    intents <- lapply(intent_paths, artifactStoreReadIntent, store = first$store)
    temporary_payloads <- vapply(
        intents,
        function(intent) intent$temporary_paths$payload,
        character(1)
    )
    expect_length(unique(temporary_payloads), 2L)
    for (intent in intents) {
        expect_identical(
            dirname(intent$temporary_paths$payload_directory),
            dirname(intent$managed_paths$payload_directory)
        )
        expect_identical(
            dirname(intent$temporary_paths$sidecar),
            dirname(intent$managed_paths$sidecar)
        )
    }

    occupied <- intents[[1L]]$managed_paths$sidecar
    occupied_path <- artifactStoreResolveFile(first$store, occupied)
    artifactStoreEnsureDirectory(first$store, dirname(occupied))
    writeLines("do not replace", occupied_path)
    expect_error(
        artifactStorePublishFile(
            first$store,
            intents[[1L]]$temporary_paths$sidecar,
            occupied
        ),
        class = "multischolar_artifact_already_exists"
    )
    expect_identical(readLines(occupied_path), "do not replace")
})

test_that("corruption and future sidecars fail read-only", {
    project_root <- withr::local_tempdir()
    fixture <- makeArtifactStoreFixture(project_root)
    ref <- artifactStoreWriteParquet(
        fixture$store,
        fixture$encoded,
        fixture$logical_key
    )
    managed <- artifactStoreManagedPaths(
        fixture$store,
        fixture$logical_key,
        ref$artifact_id
    )
    payload_path <- artifactStoreResolveFile(fixture$store, managed$payload, TRUE)
    writeBin(charToRaw("truncated"), payload_path)
    corrupt <- reconcileArtifactStore(fixture$store, repair = TRUE)
    expect_identical(corrupt$state, "corrupt_sidecar_or_payload")

    future_root <- withr::local_tempdir()
    future_fixture <- makeArtifactStoreFixture(future_root)
    future_ref <- artifactStoreWriteParquet(
        future_fixture$store,
        future_fixture$encoded,
        future_fixture$logical_key
    )
    future_paths <- artifactStoreManagedPaths(
        future_fixture$store,
        future_fixture$logical_key,
        future_ref$artifact_id
    )
    sidecar_path <- artifactStoreResolveFile(
        future_fixture$store,
        future_paths$sidecar,
        TRUE
    )
    sidecar <- jsonlite::fromJSON(
        sidecar_path,
        simplifyVector = TRUE,
        simplifyDataFrame = FALSE,
        simplifyMatrix = FALSE
    )
    sidecar$schema_version <- 2L
    jsonlite::write_json(sidecar, sidecar_path, auto_unbox = TRUE, null = "null")
    payload_digest <- artifactByteDigest(artifactStoreResolveFile(
        future_fixture$store,
        future_paths$payload,
        TRUE
    ))
    expect_error(
        artifactStoreReadSidecar(future_fixture$store, future_paths$sidecar),
        class = "multischolar_unsupported_artifact_sidecar_version"
    )
    future <- reconcileArtifactStore(future_fixture$store, repair = TRUE)
    expect_identical(future$state, "future_schema")
    expect_identical(
        artifactByteDigest(artifactStoreResolveFile(
            future_fixture$store,
            future_paths$payload,
            TRUE
        )),
        payload_digest
    )

    truncated_root <- withr::local_tempdir()
    truncated_fixture <- makeArtifactStoreFixture(truncated_root)
    truncated_ref <- artifactStoreWriteParquet(
        truncated_fixture$store,
        truncated_fixture$encoded,
        truncated_fixture$logical_key
    )
    truncated_paths <- artifactStoreManagedPaths(
        truncated_fixture$store,
        truncated_fixture$logical_key,
        truncated_ref$artifact_id
    )
    truncated_sidecar <- artifactStoreResolveFile(
        truncated_fixture$store,
        truncated_paths$sidecar,
        TRUE
    )
    writeBin(charToRaw("{\"schema\":"), truncated_sidecar)
    truncated <- reconcileArtifactStore(truncated_fixture$store, repair = TRUE)
    expect_identical(truncated$state, "corrupt_sidecar_or_payload")
})

test_that("future write intents remain read-only", {
    project_root <- withr::local_tempdir()
    fixture <- makeArtifactStoreFixture(project_root)
    expect_error(
        artifactStoreWriteParquet(
            fixture$store,
            fixture$encoded,
            fixture$logical_key,
            failure_injector = failArtifactStoreAt("before_write")
        ),
        class = "multischolar_test_artifact_failure"
    )
    intent_path <- artifactStoreIntentPaths(fixture$store)[[1L]]
    intent_file <- artifactStoreResolveFile(fixture$store, intent_path, TRUE)
    intent <- jsonlite::fromJSON(
        intent_file,
        simplifyVector = TRUE,
        simplifyDataFrame = FALSE,
        simplifyMatrix = FALSE
    )
    intent$schema_version <- 2L
    jsonlite::write_json(intent, intent_file, auto_unbox = TRUE, null = "null")
    inventory <- reconcileArtifactStore(fixture$store, repair = TRUE)
    expect_identical(inventory$state, "future_schema")
    expect_false(dir.exists(artifactPath(fixture$paths, "tables")))
})

test_that("reconciliation ignores unmanaged files and reports registry disagreement", {
    project_root <- withr::local_tempdir()
    fixture <- makeArtifactStoreFixture(project_root)
    unmanaged_root <- artifactPath(fixture$paths, "tables")
    dir.create(unmanaged_root, recursive = TRUE)
    writeLines("legacy output", file.path(unmanaged_root, "legacy-output.parquet"))
    lock_root <- artifactPath(fixture$paths, "locks")
    dir.create(lock_root, recursive = TRUE)
    stale_lock <- file.path(lock_root, "stale.lock")
    writeLines("unowned lock", stale_lock)
    expect_equal(nrow(reconcileArtifactStore(fixture$store)), 0L)
    expect_identical(readLines(stale_lock), "unowned lock")

    ref <- artifactStoreWriteParquet(
        fixture$store,
        fixture$encoded,
        fixture$logical_key
    )
    missing_id <- paste0("art_", strrep("a", 64L))
    disagreement <- reconcileArtifactStore(
        fixture$store,
        registered_artifact_ids = missing_id
    )
    expect_setequal(
        disagreement$state,
        c("committed", "registry_missing_files")
    )
    committed <- disagreement[disagreement$artifact_id == ref$artifact_id, ]
    expect_identical(committed$registry_state, "unregistered")
    expect_identical(readLines(stale_lock), "unowned lock")
})

test_that("path attacks, cross-project keys, and reserved names fail closed", {
    project_root <- withr::local_tempdir()
    fixture <- makeArtifactStoreFixture(project_root)
    bad_project <- fixture$logical_key
    bad_project$project_id <- "project-002"
    expect_error(
        artifactStoreWriteParquet(fixture$store, fixture$encoded, bad_project),
        class = "multischolar_cross_project_artifact"
    )
    for (bad_stage in c("../escape", "/absolute", "CON")) {
        bad_key <- fixture$logical_key
        bad_key$stage_id <- bad_stage
        expect_error(
            artifactStoreWriteParquet(fixture$store, fixture$encoded, bad_key),
            class = "multischolar_artifact_path_error",
            info = bad_stage
        )
    }

    external <- withr::local_tempdir()
    tables_parent <- dirname(artifactPath(fixture$paths, "tables"))
    dir.create(tables_parent, recursive = TRUE)
    linked <- file.symlink(external, artifactPath(fixture$paths, "tables"))
    if (isTRUE(linked)) {
        expect_error(
            artifactStoreWriteParquet(
                fixture$store,
                fixture$encoded,
                fixture$logical_key
            ),
            class = "multischolar_artifact_path_escape"
        )
        expect_equal(length(list.files(external, all.files = TRUE, no.. = TRUE)), 0L)
    }
})

test_that("trash is project-contained and refuses protected generations", {
    project_root <- withr::local_tempdir()
    fixture <- makeArtifactStoreFixture(project_root)
    ref <- artifactStoreWriteParquet(
        fixture$store,
        fixture$encoded,
        fixture$logical_key
    )
    expect_error(
        artifactStoreTrash(
            fixture$store,
            ref,
            current_generation_ids = ref$logical_key$generation_id
        ),
        class = "multischolar_protected_artifact_generation"
    )
    expect_error(
        artifactStoreTrash(
            fixture$store,
            ref,
            referenced_generation_ids = ref$logical_key$generation_id
        ),
        class = "multischolar_protected_artifact_generation"
    )
    receipt <- artifactStoreTrash(fixture$store, ref)
    expect_identical(receipt$status, "trashed")
    trash_path <- artifactStoreResolveFile(
        fixture$store,
        receipt$trash_relative_path,
        must_exist = TRUE
    )
    expect_true(file.exists(file.path(trash_path, "payload", "payload.parquet")))
    expect_true(file.exists(file.path(trash_path, "artifact.json")))
    expect_false(file.exists(artifactStoreResolveFile(
        fixture$store,
        ref$relative_path
    )))
})

test_that("store helpers are collated once and remain filesystem-only", {
    description <- read.dcf(testthat::test_path("..", "..", "DESCRIPTION"))
    collate <- strsplit(description[1L, "Collate"], "[[:space:]]+")[[1L]]
    collate <- gsub("^'|'$", "", collate[nzchar(collate)])
    expect_identical(sum(collate == "utils_artifact_store.R"), 1L)
    expect_identical(sum(collate == "utils_artifact_store_recovery.R"), 1L)
    expect_identical(sum(collate == "utils_artifact_store_inventory.R"), 1L)
    expect_lt(
        match("utils_artifact_store.R", collate),
        match("utils_artifact_store_recovery.R", collate)
    )
    expect_lt(
        match("utils_artifact_store_recovery.R", collate),
        match("utils_artifact_store_inventory.R", collate)
    )
    expect_lt(
        match("utils_artifact_store_inventory.R", collate),
        match("utils_workflow_context.R", collate)
    )
    sources <- unlist(lapply(
        c(
            "utils_artifact_store.R",
            "utils_artifact_store_recovery.R",
            "utils_artifact_store_inventory.R"
        ),
        function(file_name) readLines(testthat::test_path("..", "..", "R", file_name))
    ))
    expect_false(any(grepl("DuckDB|dbConnect|reactiveVal|observeEvent", sources)))
    expect_false(any(grepl("saveRDS|readRDS|chmod|Sys[.]chmod", sources)))
})
