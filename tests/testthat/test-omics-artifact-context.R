artifactContextInventory <- function() {
    jsonlite::read_json(
        testthat::test_path("..", "testdata", "omics-capabilities.json"),
        simplifyVector = FALSE
    )
}

artifactContextCapabilities <- function(inventory = artifactContextInventory()) {
    do.call(c, lapply(inventory$formats, `[[`, "capabilities"))
}

artifactContextIdentity <- function(capability, project_id = "project-1") {
    c(
        list(
            project_id = project_id,
            omic_label = paste0(capability$identity$omic_type, "_study")
        ),
        capability$identity
    )
}

artifactContextRoot <- function(prefix) {
    root <- tempfile(prefix)
    dir.create(root, recursive = TRUE)
    withr::defer(unlink(root, recursive = TRUE, force = TRUE),
        envir = parent.frame()
    )
    root
}

test_that("runtime capabilities exactly match the committed inventory", {
    expected <- artifactContextCapabilities()
    actual <- workflowCapabilityCatalogue()

    expect_identical(
        unname(vapply(actual, `[[`, character(1), "capability_id")),
        vapply(expected, `[[`, character(1), "capability_id")
    )
    expect_identical(
        unname(lapply(actual, `[[`, "identity")),
        lapply(expected, `[[`, "identity")
    )
    expect_true(all(vapply(actual, \(capability) {
        !capability$artifact_eligible &&
            !capability$auto_eligible &&
            is.null(capability$maximum_artifact_rollout)
    }, logical(1))))
    expect_identical(
        anyDuplicated(vapply(actual, \(capability) {
            workflowCapabilityKey(capability$identity)
        }, character(1))),
        0L
    )
})

test_that("all capability, backend, and project-state combinations fail closed", {
    states <- c(
        "new",
        "legacy_memory",
        "artifact_valid",
        "artifact_corrupt",
        "artifact_future_schema"
    )
    capabilities <- workflowCapabilityCatalogue()
    combinations <- 0L

    for (capability in capabilities) {
        identity <- artifactContextIdentity(capability)
        for (state in states) {
            memory <- resolveWorkflowBackend(
                identity,
                requested_backend = "memory",
                project_state = state
            )
            expect_identical(memory$effective_backend, "memory")
            expect_identical(memory$effective_rollout, "none")
            expect_identical(memory$capability_id, capability$capability_id)

            if (identical(state, "new")) {
                automatic <- resolveWorkflowBackend(
                    identity,
                    requested_backend = "auto",
                    project_state = state
                )
                expect_identical(automatic$effective_backend, "memory")
                expect_identical(
                    automatic$reason_code,
                    "auto_capability_not_promoted"
                )
            } else if (identical(state, "legacy_memory")) {
                automatic <- resolveWorkflowBackend(
                    identity,
                    requested_backend = "auto",
                    project_state = state
                )
                expect_identical(automatic$effective_backend, "memory")
                expect_identical(
                    automatic$reason_code,
                    "auto_preserve_legacy_memory"
                )
            } else {
                expected_class <- switch(state,
                    artifact_valid = "multischolar_artifact_capability_unavailable",
                    artifact_corrupt = "multischolar_corrupt_artifact_project",
                    artifact_future_schema = "multischolar_future_artifact_schema"
                )
                expect_error(
                    resolveWorkflowBackend(
                        identity,
                        requested_backend = "auto",
                        project_state = state
                    ),
                    class = expected_class
                )
            }
            expect_error(
                resolveWorkflowBackend(
                    identity,
                    requested_backend = "artifact",
                    project_state = state
                ),
                class = "multischolar_artifact_not_certified"
            )
            combinations <- combinations + 3L
        }
    }
    expect_identical(combinations, length(capabilities) * length(states) * 3L)
})

test_that("promoted capabilities obey rollout and migration bounds", {
    promoted <- workflowCapabilityCatalogue()[[1L]]
    promoted$artifact_eligible <- TRUE
    promoted$auto_eligible <- TRUE
    promoted$maximum_artifact_rollout <- "read_through"
    capabilities <- list(promoted)
    identity <- artifactContextIdentity(promoted)

    automatic <- resolveWorkflowBackend(
        identity,
        requested_backend = "auto",
        requested_rollout = "evict",
        project_state = "new",
        capabilities = capabilities
    )
    expect_named(automatic, c(
        "requested_backend", "effective_backend", "effective_rollout",
        "capability_id", "capability_version", "reason_code", "project_state"
    ))
    expect_identical(automatic$effective_backend, "artifact")
    expect_identical(automatic$effective_rollout, "read_through")
    expect_identical(automatic$reason_code, "auto_promoted_new_project")

    legacy <- resolveWorkflowBackend(
        identity,
        requested_backend = "auto",
        project_state = "legacy_memory",
        capabilities = capabilities
    )
    expect_identical(legacy$effective_backend, "memory")
    expect_identical(legacy$reason_code, "auto_preserve_legacy_memory")

    migrated <- resolveWorkflowBackend(
        identity,
        requested_backend = "auto",
        project_state = "legacy_memory",
        migration_requested = TRUE,
        capabilities = capabilities
    )
    expect_identical(migrated$effective_backend, "artifact")
    expect_identical(migrated$reason_code, "auto_migration_requested")

    existing <- resolveWorkflowBackend(
        identity,
        requested_backend = "auto",
        project_state = "artifact_valid",
        capabilities = capabilities
    )
    expect_identical(existing$effective_backend, "artifact")
    expect_identical(existing$reason_code, "auto_existing_artifact")

    expect_error(
        resolveWorkflowBackend(
            identity,
            requested_backend = "artifact",
            requested_rollout = "evict",
            capabilities = capabilities
        ),
        class = "multischolar_artifact_rollout_not_certified"
    )
    expect_error(
        resolveWorkflowBackend(
            identity,
            requested_backend = "artifact",
            project_state = "legacy_memory",
            capabilities = capabilities
        ),
        class = "multischolar_artifact_migration_required"
    )

    unknown <- identity
    unknown$acquisition_mode <- "unknown"
    expect_identical(
        resolveWorkflowBackend(unknown, requested_backend = "auto")$effective_backend,
        "memory"
    )
    expect_error(
        resolveWorkflowBackend(unknown, requested_backend = "artifact"),
        class = "multischolar_unknown_workflow_capability"
    )
})

test_that("unsupported forced artifact mode fails before project resources", {
    root <- artifactContextRoot("artifact-context-preflight-")
    calls <- new.env(parent = emptyenv())
    calls$path_builder <- 0L
    calls$state_detector <- 0L
    context <- createWorkflowContext(
        experiment_paths = list(base_dir = root),
        omic_type = "proteomics",
        storage_policy = list(requested_backend = "artifact")
    )

    expect_error(
        bindWorkflowContextFromImport(
            context,
            workflow_type = "DIA",
            input_format = "diann",
            data_level = "peptide",
            path_builder = function(...) {
                calls$path_builder <- calls$path_builder + 1L
                stop("path builder must not run")
            },
            project_state_detector = function(...) {
                calls$state_detector <- calls$state_detector + 1L
                stop("state detector must not run")
            }
        ),
        class = "multischolar_artifact_not_certified"
    )
    expect_identical(calls$path_builder, 0L)
    expect_identical(calls$state_detector, 0L)
    expect_false(dir.exists(file.path(root, "state")))
})

test_that("workflow context binds once and notifies early and late consumers", {
    root <- artifactContextRoot("artifact-context-memory-")
    context <- createWorkflowContext(
        experiment_paths = list(base_dir = root, omic_label = "proteomics_study"),
        omic_type = "proteomics",
        storage_policy = list(requested_backend = "memory")
    )
    early <- list()
    unsubscribe <- context$observeBinding(function(snapshot) {
        early[[length(early) + 1L]] <<- snapshot
    })

    snapshot <- bindWorkflowContextFromImport(
        context,
        workflow_type = "DIA",
        input_format = "diann",
        data_level = "peptide"
    )
    expect_true(context$isBound())
    expect_length(early, 1L)
    expect_identical(early[[1L]], snapshot)
    expect_identical(snapshot$resolution$effective_backend, "memory")
    expect_null(snapshot$path_metadata)
    expect_false(dir.exists(file.path(root, "state")))

    late <- NULL
    context$observeBinding(function(value) late <<- value)
    expect_identical(late, snapshot)
    expect_false(context$bind(
        context$getIdentity(),
        context$getStorageDecision(),
        context$getPaths()
    ))
    changed <- context$getIdentity()
    changed$input_format <- "changed"
    expect_error(
        context$bind(changed, context$getStorageDecision()),
        class = "multischolar_workflow_context_already_bound"
    )
    unsubscribe()
})

test_that("explicit memory binding does not inspect artifact resources", {
    missing_root <- tempfile("artifact-context-memory-missing-")
    context <- createWorkflowContext(
        experiment_paths = list(base_dir = missing_root),
        omic_type = "proteomics",
        storage_policy = list(requested_backend = "memory")
    )
    detector_calls <- 0L

    bindWorkflowContextFromImport(
        context,
        workflow_type = "DIA",
        input_format = "diann",
        data_level = "peptide",
        project_state_detector = function(...) {
            detector_calls <<- detector_calls + 1L
            stop("memory mode must not inspect artifact resources")
        }
    )
    expect_identical(detector_calls, 0L)
    expect_identical(
        context$getStorageDecision()$effective_backend,
        "memory"
    )
    expect_false(dir.exists(missing_root))
})

test_that("legacy workflow path lists infer one project root without mutation", {
    root <- artifactContextRoot("artifact-context-legacy-paths-")
    experiment_paths <- list(
        data_dir = file.path(root, "data", "proteomics"),
        results_dir = file.path(root, "results", "proteomics_study"),
        source_dir = file.path(root, "scripts", "proteomics_study")
    )
    context <- createWorkflowContext(
        experiment_paths = experiment_paths,
        omic_type = "proteomics",
        experiment_label = "proteomics_study"
    )

    expect_identical(context$getProjectRoot(), normalizePath(root))
    expect_false(context$isBound())
    expect_false(dir.exists(file.path(root, "state")))
    expect_identical(names(experiment_paths), c(
        "data_dir", "results_dir", "source_dir"
    ))
})

test_that("artifact paths are deterministic, relative, labelled, and contained", {
    root <- artifactContextRoot("artifact-paths-")
    capability <- workflowCapabilityCatalogue()[[1L]]
    identity <- artifactContextIdentity(capability)
    paths <- buildArtifactPaths(root, identity)
    repeated <- buildArtifactPaths(root, identity)
    metadata <- artifactPathMetadata(paths)

    expect_identical(metadata, artifactPathMetadata(repeated))
    expect_identical(metadata$labels$omic_label, "proteomics_study")
    expect_identical(metadata$labels$workflow_slug, "prot_dia")
    expect_identical(
        metadata$relative_paths$tables,
        "data/proteomics_study/prot_dia/tables"
    )
    expect_identical(
        metadata$legacy_read_fallbacks$data_root$relative_path,
        "data/proteomics"
    )
    expect_identical(
        metadata$legacy_read_fallbacks$data_root$access,
        "read_only"
    )
    expect_false(any(startsWith(unlist(metadata$relative_paths), "/")))
    expect_identical(
        artifactPath(paths, "tables"),
        file.path(normalizePath(root), "data", "proteomics_study", "prot_dia", "tables")
    )
    expect_false(dir.exists(file.path(root, "data")))

    expect_error(
        artifactResolveContainedPath(root, "../escape"),
        class = "multischolar_invalid_relative_artifact_path"
    )
    expect_error(
        artifactResolveContainedPath(root, normalizePath(root)),
        class = "multischolar_invalid_relative_artifact_path"
    )
    invalid_identity <- identity
    invalid_identity$omic_label <- "CON"
    expect_error(
        buildArtifactPaths(root, invalid_identity),
        class = "multischolar_invalid_artifact_path_component"
    )

    outside <- artifactContextRoot("artifact-paths-outside-")
    dir.create(file.path(root, "data"), recursive = TRUE)
    linked <- file.symlink(outside, file.path(root, "data", "escaped"))
    if (isTRUE(linked)) {
        expect_error(
            artifactResolveContainedPath(root, "data/escaped/workflow/tables"),
            class = "multischolar_artifact_path_escape"
        )
    }
})

test_that("project-state detection distinguishes legacy and artifact evidence", {
    identity <- artifactContextIdentity(workflowCapabilityCatalogue()[[1L]])

    new_root <- artifactContextRoot("artifact-state-new-")
    expect_identical(detectWorkflowProjectState(new_root, identity)$status, "new")

    legacy_root <- artifactContextRoot("artifact-state-legacy-")
    dir.create(file.path(legacy_root, "data", "proteomics"), recursive = TRUE)
    writeLines("legacy", file.path(legacy_root, "data", "proteomics", "input.tsv"))
    expect_identical(
        detectWorkflowProjectState(legacy_root, identity)$status,
        "legacy_memory"
    )

    artifact_root <- artifactContextRoot("artifact-state-valid-")
    artifact_paths <- buildArtifactPaths(artifact_root, identity)
    manifest <- artifactPath(artifact_paths, "artifact_manifest")
    dir.create(dirname(manifest), recursive = TRUE)
    jsonlite::write_json(list(
        schema = "multischolar.artifact_manifest",
        schema_version = 1L
    ), manifest, auto_unbox = TRUE)
    expect_identical(
        detectWorkflowProjectState(artifact_root, identity)$status,
        "artifact_valid"
    )

    future_root <- artifactContextRoot("artifact-state-future-")
    future_paths <- buildArtifactPaths(future_root, identity)
    future_manifest <- artifactPath(future_paths, "artifact_manifest")
    dir.create(dirname(future_manifest), recursive = TRUE)
    jsonlite::write_json(list(
        schema = "multischolar.artifact_manifest",
        schema_version = 2L
    ), future_manifest, auto_unbox = TRUE)
    expect_identical(
        detectWorkflowProjectState(future_root, identity)$status,
        "artifact_future_schema"
    )

    corrupt_root <- artifactContextRoot("artifact-state-corrupt-")
    corrupt_paths <- buildArtifactPaths(corrupt_root, identity)
    corrupt_manifest <- artifactPath(corrupt_paths, "artifact_manifest")
    dir.create(dirname(corrupt_manifest), recursive = TRUE)
    writeLines("{broken", corrupt_manifest)
    expect_identical(
        detectWorkflowProjectState(corrupt_root, identity)$status,
        "artifact_corrupt"
    )

    evidence_root <- artifactContextRoot("artifact-state-evidence-")
    evidence_paths <- buildArtifactPaths(evidence_root, identity)
    tables <- artifactPath(evidence_paths, "tables")
    dir.create(tables, recursive = TRUE)
    writeLines("orphan", file.path(tables, "table.parquet"))
    expect_identical(
        detectWorkflowProjectState(evidence_root, identity)$status,
        "artifact_corrupt"
    )

    shared_root <- artifactContextRoot("artifact-state-shared-registry-")
    dir.create(file.path(shared_root, "state"), recursive = TRUE)
    file.create(file.path(shared_root, "state", "multischolar.duckdb"))
    expect_identical(
        detectWorkflowProjectState(shared_root, identity)$status,
        "new"
    )
})

test_that("artifact helpers are collated once before every caller", {
    description <- read.dcf(
        testthat::test_path("..", "..", "DESCRIPTION"),
        fields = "Collate"
    )[[1L]]
    collate <- strsplit(description, "[[:space:]]+")[[1L]]
    collate <- gsub("^'|'$", "", collate[nzchar(collate)])
    helpers <- c(
        "utils_workflow_capabilities.R",
        "utils_artifact_paths.R",
        "utils_workflow_context.R"
    )
    callers <- c(
        "func_general_filemgmt.R",
        "mod_proteomics.R",
        "mod_metabolomics.R",
        "mod_lipidomics.R"
    )

    expect_true(all(vapply(helpers, \(helper) sum(collate == helper) == 1L, logical(1))))
    expect_true(max(match(helpers, collate)) < min(match(callers, collate)))
})
