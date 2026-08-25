autoPolicyReadRecord <- function() {
    jsonlite::read_json(
        operationalArtifactRepoPath(
            "tests", "testdata", "omics-parity",
            "all-omics-auto-policy-v1.json"
        ),
        simplifyVector = FALSE
    )
}

autoPolicyIdentity <- function(capability) {
    c(
        list(
            project_id = "auto-policy-project",
            omic_label = paste0(capability$identity$omic_type, "-study")
        ),
        capability$identity
    )
}

autoPolicyDescriptionPackages <- function(field) {
    description <- read.dcf(operationalArtifactRepoPath("DESCRIPTION"))
    values <- strsplit(description[[1L, field]], ",", fixed = TRUE)[[1L]]
    trimws(sub("[[:space:]]*[(].*$", "", values))
}

test_that("auto policy is bound to the immutable empty closeout set", {
    policy <- autoPolicyReadRecord()
    closeout_path <- operationalArtifactRepoPath(policy$closeout_path)
    expect_identical(
        policy$schema,
        "multischolar.all_omics_auto_policy"
    )
    expect_identical(policy$schema_version, "1.0.0")
    expect_identical(policy$decision, "no_promotion")
    expect_true(file.exists(closeout_path))
    expect_identical(
        artifactByteDigest(closeout_path),
        policy$closeout_sha256
    )
    closeout <- operationalArtifactReadCloseout()
    expect_length(closeout$promoted_capability_ids, 0L)
    expect_length(policy$promoted_capability_ids, 0L)
    expect_identical(policy$default_requested_backend, "auto")
    expect_identical(
        policy$effective_backend_for_new_projects,
        "memory"
    )
    expect_identical(
        normalizeWorkflowStoragePolicy()$requested_backend,
        "auto"
    )
})

test_that("every exact capability resolves auto and memory fail closed", {
    capabilities <- mergeWorkflowDescriptorCapabilities()
    descriptors <- artifactWorkflowDescriptorCatalogue()$descriptors
    policy <- autoPolicyReadRecord()
    canary_ids <- unlist(
        policy$forced_canary_capability_ids,
        use.names = FALSE
    )
    expect_setequal(canary_ids, names(descriptors))

    for (capability in capabilities) {
        identity <- autoPolicyIdentity(capability)
        automatic <- resolveWorkflowBackend(
            identity,
            requested_backend = "auto",
            project_state = "new",
            capabilities = capabilities
        )
        expect_identical(automatic$effective_backend, "memory")
        expect_identical(
            automatic$reason_code,
            "auto_capability_not_promoted"
        )
        memory <- resolveWorkflowBackend(
            identity,
            requested_backend = "memory",
            project_state = "new",
            capabilities = capabilities
        )
        expect_identical(memory$effective_backend, "memory")
        legacy <- resolveWorkflowBackend(
            identity,
            requested_backend = "auto",
            project_state = "legacy_memory",
            capabilities = capabilities
        )
        expect_identical(legacy$effective_backend, "memory")
        expect_identical(legacy$reason_code, "auto_preserve_legacy_memory")

        if (capability$capability_id %in% canary_ids) {
            forced <- resolveWorkflowBackend(
                identity,
                requested_backend = "artifact",
                requested_rollout = "dual_write",
                project_state = "new",
                migration_requested = TRUE,
                capabilities = capabilities
            )
            expect_identical(forced$effective_backend, "artifact")
            expect_identical(forced$effective_rollout, "dual_write")
            existing <- resolveWorkflowBackend(
                identity,
                requested_backend = "auto",
                project_state = "artifact_valid",
                capabilities = capabilities
            )
            expect_identical(existing$effective_backend, "artifact")
            expect_identical(existing$reason_code, "auto_existing_artifact")
        } else {
            expect_error(
                resolveWorkflowBackend(
                    identity,
                    requested_backend = "artifact",
                    project_state = "new",
                    migration_requested = TRUE,
                    capabilities = capabilities
                ),
                class = "multischolar_artifact_not_certified"
            )
            expect_error(
                resolveWorkflowBackend(
                    identity,
                    requested_backend = "auto",
                    project_state = "artifact_valid",
                    capabilities = capabilities
                ),
                class = "multischolar_artifact_capability_unavailable"
            )
        }
        expect_error(
            resolveWorkflowBackend(
                identity,
                requested_backend = "auto",
                project_state = "artifact_corrupt",
                capabilities = capabilities
            ),
            class = "multischolar_corrupt_artifact_project"
        )
        expect_error(
            resolveWorkflowBackend(
                identity,
                requested_backend = "auto",
                project_state = "artifact_future_schema",
                capabilities = capabilities
            ),
            class = "multischolar_future_artifact_schema"
        )
    }

    unknown <- autoPolicyIdentity(capabilities[[1L]])
    unknown$input_format <- "unknown"
    expect_identical(
        resolveWorkflowBackend(
            unknown,
            requested_backend = "auto",
            capabilities = capabilities
        )$effective_backend,
        "memory"
    )
    expect_error(
        resolveWorkflowBackend(
            unknown,
            requested_backend = "artifact",
            capabilities = capabilities
        ),
        class = "multischolar_unknown_workflow_capability"
    )
})

test_that("optional artifact dependency policy is deliberate and actionable", {
    policy <- autoPolicyReadRecord()
    imports <- autoPolicyDescriptionPackages("Imports")
    suggests <- autoPolicyDescriptionPackages("Suggests")
    dependency_policy <- policy$dependency_policy
    dependencies <- vapply(
        dependency_policy,
        `[[`,
        character(1),
        "package"
    )
    expect_setequal(dependencies, c("arrow", "DBI", "duckdb", "filelock"))
    expect_true("arrow" %in% imports)
    expect_false(any(c("DBI", "duckdb", "filelock") %in% imports))
    expect_true(all(c("DBI", "duckdb", "filelock") %in% suggests))
    expect_setequal(
        PROJECT_REGISTRY_REQUIRED_PACKAGES,
        c("DBI", "duckdb", "filelock")
    )
    for (package in PROJECT_REGISTRY_REQUIRED_PACKAGES) {
        availability <- stats::setNames(
            rep(TRUE, length(PROJECT_REGISTRY_REQUIRED_PACKAGES)),
            PROJECT_REGISTRY_REQUIRED_PACKAGES
        )
        availability[[package]] <- FALSE
        error <- rlang::catch_cnd(
            projectRegistryAssertDependencies(availability)
        )
        expect_s3_class(error, "multischolar_missing_artifact_dependencies")
        expect_identical(error$missing_packages, package)
        expect_match(conditionMessage(error), package, fixed = TRUE)
    }
    namespace <- readLines(operationalArtifactRepoPath("NAMESPACE"), warn = FALSE)
    expect_false(any(namespace %in% c(
        "import(DBI)", "import(duckdb)", "import(filelock)"
    )))
})

test_that("memory contexts remain dependency inert and unmigrated", {
    root <- withr::local_tempdir(pattern = "omics-art-049-memory-")
    context <- createWorkflowContext(
        list(base_dir = root, omic_label = "memory-study"),
        "proteomics",
        "memory-study",
        storage_policy = list(
            requested_backend = "memory",
            project_id = "auto-policy-memory"
        )
    )
    detector_calls <- 0L
    bindWorkflowContextFromImport(
        context,
        workflow_type = "DIA",
        input_format = "diann",
        data_level = "peptide",
        descriptor_catalogue = artifactWorkflowDescriptorCatalogue(),
        project_state_detector = \(...) {
            detector_calls <<- detector_calls + 1L
            testthat::fail("memory binding inspected artifact state")
        }
    )
    expect_identical(detector_calls, 0L)
    expect_identical(
        context$getStorageDecision()$effective_backend,
        "memory"
    )
    expect_null(context$getPaths())
    expect_null(projectRegistryForContext(context))
    manager <- newWorkflowState(workflow_context = context)
    expect_s3_class(manager, "WorkflowState")
    expect_false(inherits(manager, "ArtifactWorkflowState"))
})

test_that("support auxiliary and compatibility policies remain unchanged", {
    policy <- autoPolicyReadRecord()
    closeout <- operationalArtifactReadCloseout()
    capability_ids <- names(workflowCapabilityCatalogue())
    expect_setequal(
        unlist(closeout$capability_ids, use.names = FALSE),
        capability_ids
    )
    expect_false(any(grepl(
        "transcript|integration|phospho|multiomics",
        capability_ids,
        ignore.case = TRUE
    )))
    expect_true(is.function(processMultisiteEvidence))
    expect_true(is.function(plotMofaWeights))
    compatibility <- policy$compatibility_policy
    retained <- setdiff(names(compatibility), c(
        "removal_schedule_status", "earliest_removal_release"
    ))
    expect_true(all(vapply(
        compatibility[retained],
        identical,
        logical(1),
        "retained"
    )))
    expect_identical(
        compatibility$removal_schedule_status,
        "not_started_no_promoted_tuple"
    )
    expect_null(compatibility$earliest_removal_release)
    expect_identical(policy$rollback$new_work_backend, "memory")
    expect_false(policy$rollback$delete_immutable_artifacts)
})

test_that("janitor promotion bands remain closed", {
    policy <- autoPolicyReadRecord()
    runtime_files <- list.files(
        operationalArtifactRepoPath("R"),
        pattern = "[.]R$",
        full.names = TRUE
    )
    line_counts <- vapply(
        runtime_files,
        \(path) length(readLines(path, warn = FALSE)),
        integer(1)
    )
    expect_lte(
        max(line_counts),
        as.integer(policy$janitor_policy$maximum_runtime_file_lines)
    )
    expect_identical(
        unname(line_counts[basename(runtime_files) ==
            "utils_workflow_state_artifact_backend.R"]),
        1000L
    )
    audit <- paste(readLines(
        operationalArtifactRepoPath(
            "tools", "refactor", "FULL_FUNCTION_PARITY_AUDIT.md"
        ),
        warn = FALSE
    ), collapse = "\n")
    expect_match(audit, "zero duplicate entity keys", fixed = TRUE)
    expect_match(audit, "zero stale extraction", fixed = TRUE)
    expect_match(audit, "zero filename violations", fixed = TRUE)
})
