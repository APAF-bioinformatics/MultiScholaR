hybridLoadAssignments <- function(path, names) {
    expressions <- parse(testthat::test_path(path))
    target <- parent.frame()
    for (expression in expressions) {
        assignment <- is.call(expression) &&
            identical(expression[[1L]], as.name("<-"))
        if (isTRUE(assignment) && as.character(expression[[2L]]) %in% names) {
            eval(expression, envir = target)
        }
    }
    invisible(TRUE)
}

hybridLoadAssignments(
    "test-omics-artifact-dia-normalization.R",
    c(
        "omics014Capability", "omics014Context", "omics014Manager",
        "omics014Protein"
    )
)
hybridLoadAssignments(
    "test-omics-artifact-metab-closeout.R",
    c("metab037Context", "metab037Manager")
)
hybridLoadAssignments(
    "test-omics-artifact-lipid-closeout.R",
    c("lipid046Context", "lipid046Manager")
)
hybridLoadAssignments(
    "test-omics-artifact-prot-nondia-norm.R",
    c(
        "nondiaArtifact025RepoPath", "nondiaArtifact025Scenarios",
        "nondiaArtifact025Importer", "nondiaArtifact025WorkflowType",
        "nondiaArtifact025VendorObject", "nondiaArtifact025CapabilityId",
        "nondiaArtifact025Manager"
    )
)

hybridMemoryContext <- function(root, project_id) {
    capability <- omics014Capability()
    capability$artifact_eligible <- FALSE
    capability$auto_eligible <- FALSE
    context <- createWorkflowContext(
        list(
            base_dir = root,
            project_id = project_id,
            omic_label = "memory-proteomics"
        ),
        "proteomics",
        "memory-proteomics",
        storage_policy = list(
            requested_backend = "memory",
            requested_rollout = "dual_write",
            migration_requested = FALSE,
            project_id = project_id
        )
    )
    bindWorkflowContextFromImport(
        context,
        workflow_type = "DIA",
        input_format = "diann",
        data_level = "peptide",
        capabilities = list(capability)
    )
    context
}

hybridNondiaContext <- function(root, descriptor, project_id, backend) {
    context <- createWorkflowContext(
        list(
            base_dir = root,
            project_id = project_id,
            omic_type = "proteomics",
            omic_label = "nondia-study"
        ),
        "proteomics",
        "nondia-study",
        storage_policy = list(
            requested_backend = backend,
            requested_rollout = "dual_write",
            migration_requested = identical(backend, "artifact"),
            project_id = project_id
        )
    )
    identity <- descriptor$identity
    bindWorkflowContextFromImport(
        context,
        workflow_type = identity$workflow_type,
        input_format = identity$input_format,
        data_level = identity$data_level,
        descriptor_catalogue = artifactWorkflowDescriptorCatalogue()
    )
    context
}

hybridMemoryManager <- function(context, object, state_name, workflow_type) {
    manager <- newWorkflowState(workflow_context = context)
    manager$setWorkflowType(workflow_type)
    manager$saveState(
        state_name,
        object,
        object@args,
        "hybrid memory"
    )
    manager
}

hybridExpectClosed <- function(manager) {
    info <- manager$getResourceInfo()
    expect_true(info$closed)
    expect_false(info$registry_connection)
    expect_false(info$writer_guard)
    expect_identical(info$hydration_cache_entries, 0L)
}

test_that("real hybrid coordinators isolate roots writers state and rollback", {
    for (package in c("DBI", "duckdb", "filelock", "arrow")) {
        skip_if_not_installed(package)
    }
    root <- withr::local_tempdir(pattern = "omics-art-047-hybrid-")
    globals <- c(
        "data_tbl", "data_cln", "design_matrix", "contrasts_tbl",
        "uniprot_dat_cln", "filtering_progress", "project_dirs",
        "config_list", "filtering_progress_metabolomics",
        "filtering_progress_lipidomics", "workflow_context"
    )
    globals_before <- vapply(
        globals,
        exists,
        logical(1),
        envir = .GlobalEnv,
        inherits = FALSE
    )

    project_id <- "hybrid-project"
    nondia_scenario <- Filter(
        \(scenario) identical(scenario$format, "maxquant"),
        nondiaArtifact025Scenarios()
    )[[1L]]
    nondia_descriptor <- artifactProteomicsNonDiaWorkflowDescriptors()[[
        nondiaArtifact025CapabilityId(nondia_scenario$format)
    ]]
    dia_context <- omics014Context(root, project_id)
    nondia_context <- hybridNondiaContext(
        root,
        nondia_descriptor,
        project_id,
        "artifact"
    )
    metab_context <- metab037Context(root, project_id)
    lipid_context <- lipid046Context(root, project_id)
    dia_memory_context <- hybridMemoryContext(root, project_id)
    nondia_memory_context <- hybridNondiaContext(
        root,
        nondia_descriptor,
        project_id,
        "memory"
    )
    dia_object <- omics014Protein()
    nondia_object <- suppressMessages(
        nondiaArtifact025VendorObject(nondia_scenario)
    )
    metab_object <- module_ci_metab_norm_object(
        layout = "combined",
        positive_only = TRUE
    )
    lipid_object <- module_ci_lipid_norm_object(
        layout = "combined",
        positive_only = TRUE
    )

    dia_memory_manager <- hybridMemoryManager(
        dia_memory_context,
        dia_object,
        "memory_protein",
        "DIA"
    )
    nondia_memory_manager <- hybridMemoryManager(
        nondia_memory_context,
        nondia_object,
        "memory_nondia_protein",
        nondia_descriptor$identity$workflow_type
    )
    dia_memory_snapshot <- dia_memory_manager$exportState()
    nondia_memory_snapshot <- nondia_memory_manager$exportState()

    contexts <- list(
        dia = dia_context,
        nondia = nondia_context,
        metabolomics = metab_context,
        lipidomics = lipid_context,
        dia_memory = dia_memory_context,
        nondia_memory = nondia_memory_context
    )
    roots <- vapply(
        contexts,
        \(context) context$getProjectRoot(),
        character(1)
    )
    expect_true(all(roots == roots[[1L]]))
    workflow_ids <- vapply(
        contexts,
        \(context) context$getIdentity()$workflow_id,
        character(1)
    )
    expect_identical(workflow_ids[["dia"]], workflow_ids[["nondia"]])
    expect_identical(workflow_ids[["dia"]], workflow_ids[["dia_memory"]])
    expect_identical(
        workflow_ids[["nondia"]],
        workflow_ids[["nondia_memory"]]
    )
    expect_true(all(vapply(
        contexts,
        \(context) identical(context$getIdentity()$project_id, project_id),
        logical(1)
    )))
    artifact_names <- c("dia", "nondia", "metabolomics", "lipidomics")
    memory_names <- c("dia_memory", "nondia_memory")
    capability_keys <- vapply(
        contexts[artifact_names],
        \(context) workflowCapabilityKey(context$getIdentity()),
        character(1)
    )
    expect_identical(anyDuplicated(capability_keys), 0L)
    expect_true(all(vapply(
        contexts[artifact_names],
        \(context) identical(
            context$getStorageDecision()$effective_backend,
            "artifact"
        ),
        logical(1)
    )))
    expect_true(all(vapply(
        contexts[memory_names],
        \(context) identical(
            context$getStorageDecision()$effective_backend,
            "memory"
        ) && is.null(context$getPaths()),
        logical(1)
    )))

    isolated_paths <- c(
        "data_root", "cache", "staging", "trash",
        "workflow_state_root", "generations", "artifact_manifest"
    )
    for (path_name in isolated_paths) {
        relative_paths <- vapply(
            contexts[artifact_names],
            \(context) context$getPaths()$relative_paths[[path_name]],
            character(1)
        )
        expect_identical(anyDuplicated(relative_paths), 0L)
    }
    expected_lock_path <- file.path(
        artifactPath(dia_context$getPaths(), "locks"),
        "project-registry.lock"
    )
    expect_false(file.exists(expected_lock_path))

    dia_manager <- omics014Manager(dia_context)
    dia_manager$saveState(
        "protein_replicate_filtered",
        dia_object,
        dia_object@args,
        "hybrid DIA"
    )
    expect_identical(workflowStateType(dia_manager), "DIA")
    expect_s3_class(dia_memory_manager, "WorkflowState")
    expect_s3_class(nondia_memory_manager, "WorkflowState")
    expect_false(inherits(dia_memory_manager, "ArtifactWorkflowState"))
    expect_false(inherits(nondia_memory_manager, "ArtifactWorkflowState"))

    lock_paths <- vapply(
        contexts[artifact_names],
        \(context) file.path(
            artifactPath(context$getPaths(), "locks"),
            "project-registry.lock"
        ),
        character(1)
    )
    expect_identical(length(unique(lock_paths)), 1L)
    expect_true(file.exists(lock_paths[[1L]]))
    expect_error(
        omics014Manager(dia_context),
        class = "multischolar_registry_writer_busy"
    )
    expect_error(
        metab037Manager(metab_context, metab_object),
        class = "multischolar_registry_writer_busy"
    )
    expect_identical(dia_memory_manager$getState(), dia_object)
    expect_identical(nondia_memory_manager$getState(), nondia_object)
    expect_identical(dia_manager$getCacheInfo()$entries, 1L)
    expect_true(dia_manager$releaseCache())
    expect_identical(dia_manager$getCacheInfo()$entries, 0L)
    expect_identical(dia_manager$getState(), dia_object)
    expect_identical(dia_manager$getCacheInfo()$entries, 1L)
    dia_generation <- dia_manager$getCurrentGenerationId()
    dia_snapshot <- dia_manager$exportState()
    expect_identical(dia_snapshot$workflow_id, workflow_ids[["dia"]])
    expect_identical(dia_snapshot$current_generation_id, dia_generation)
    expect_true(dia_manager$close())
    hybridExpectClosed(dia_manager)

    nondia_manager <- nondiaArtifact025Manager(
        nondia_context,
        nondia_descriptor,
        nondia_object
    )
    nondia_generation <- nondia_manager$getCurrentGenerationId()
    expect_identical(
        workflowStateType(nondia_manager),
        nondia_descriptor$identity$workflow_type
    )
    expect_identical(nondia_memory_manager$getState(), nondia_object)
    expect_true(nondia_manager$close())
    hybridExpectClosed(nondia_manager)

    metab_manager <- metab037Manager(metab_context, metab_object)
    expect_identical(workflowStateType(metab_manager), "metabolomics_standard")
    metab_generation <- metab_manager$getCurrentGenerationId()
    expect_identical(dia_memory_manager$getState(), dia_object)
    expect_true(metab_manager$close())
    hybridExpectClosed(metab_manager)

    lipid_manager <- lipid046Manager(lipid_context, lipid_object)
    expect_identical(workflowStateType(lipid_manager), "lipidomics_standard")
    lipid_generation <- lipid_manager$getCurrentGenerationId()
    expect_error(
        lipid_manager$saveState(
            "invalid_lipid_state",
            dia_object,
            list(),
            "wrong omic"
        )
    )
    expect_identical(
        lipid_manager$getCurrentGenerationId(),
        lipid_generation
    )

    lipid_manager$saveState(
        "lipid_second",
        lipid_object,
        lipid_object@args,
        "hybrid lipid second"
    )
    lipid_manager$revertToState("lipid_norm_complete")
    expect_identical(lipid_manager$getState(), lipid_object)
    expect_identical(dia_memory_manager$getState(), dia_object)
    expect_identical(nondia_memory_manager$getState(), nondia_object)
    expect_identical(dia_memory_manager$exportState(), dia_memory_snapshot)
    expect_identical(
        nondia_memory_manager$exportState(),
        nondia_memory_snapshot
    )
    expect_true(lipid_manager$close())
    hybridExpectClosed(lipid_manager)

    dia_reopened <- omics014Manager(dia_context)
    expect_identical(dia_reopened$getCurrentGenerationId(), dia_generation)
    expect_identical(dia_reopened$getState(), dia_object)
    expect_true(dia_reopened$close())
    hybridExpectClosed(dia_reopened)
    nondia_reopened <- newWorkflowState(
        workflow_context = nondia_context,
        workflow_descriptor = nondia_descriptor,
        descriptor_catalogue = artifactWorkflowDescriptorCatalogue(),
        codec_catalogue = artifactS4CodecCatalogue()
    )
    expect_identical(
        nondia_reopened$getCurrentGenerationId(),
        nondia_generation
    )
    expect_identical(nondia_reopened$getState(), nondia_object)
    expect_true(nondia_reopened$close())
    hybridExpectClosed(nondia_reopened)
    metab_reopened <- ArtifactWorkflowState$new(
        workflow_context = metab_context,
        dehydrate_fn = dehydrateMetabolomicsS4Artifact,
        validate_bundle_fn = validateMetabolomicsS4Bundle,
        hydrate_fn = hydrateMetabolomicsS4Artifact,
        descriptor_contract = artifactStageDescriptorContract(
            artifactMetabolomicsWorkflowDescriptor()
        )
    )
    expect_identical(metab_reopened$getCurrentGenerationId(), metab_generation)
    expect_identical(metab_reopened$getState(), metab_object)
    expect_true(metab_reopened$close())
    hybridExpectClosed(metab_reopened)
    expect_identical(dia_memory_manager$getState(), dia_object)
    expect_identical(nondia_memory_manager$getState(), nondia_object)

    registry <- projectRegistryForContext(dia_context)
    registry_ids <- vapply(contexts[artifact_names], \(context) {
        identity <- context$getIdentity()
        store <- newArtifactStore(context$getPaths(), identity$project_id)
        scope_path <- artifactStoreResolveFile(
            store,
            artifactRegistryScopePath(store),
            must_exist = TRUE
        )
        expect_true(file.exists(scope_path))
        artifactRegistryIdentity(store, identity)$workflow_id
    }, character(1))
    expect_identical(anyDuplicated(registry_ids), 0L)
    inspection <- openProjectRegistryReadOnly(registry)
    on.exit(closeProjectRegistry(inspection), add = TRUE)
    workflows <- projectRegistryQuery(inspection, "workflows")
    expect_setequal(workflows$workflow_id, registry_ids)
    expect_identical(anyDuplicated(workflows$workflow_id), 0L)
    expect_error(
        projectRegistryQuery(
            inspection,
            "states",
            list(workflow_id = workflow_ids[["dia"]])
        ),
        class = "multischolar_ambiguous_registry_workflow_id"
    )
    expected_states <- c(
        dia = "protein_replicate_filtered",
        nondia = "protein_s4_initial",
        metabolomics = "metab_norm_complete",
        lipidomics = "lipid_norm_complete"
    )
    expected_generations <- c(
        dia = dia_generation,
        nondia = nondia_generation,
        metabolomics = metab_generation,
        lipidomics = lipid_generation
    )
    for (name in artifact_names) {
        rows <- projectRegistryQuery(
            inspection,
            "states",
            list(workflow_id = registry_ids[[name]])
        )
        current <- rows[rows$status == "current", , drop = FALSE]
        expect_identical(nrow(current), 1L)
        expect_identical(current$logical_name, expected_states[[name]])
        expect_identical(
            current$generation_id,
            expected_generations[[name]]
        )
        prefix <- contexts[[name]]$getPaths()$relative_paths$generations
        expect_true(startsWith(current$manifest_relative_path, prefix))
    }
    lipid_rows <- projectRegistryQuery(
        inspection,
        "states",
        list(workflow_id = registry_ids[["lipidomics"]])
    )
    expect_true(any(lipid_rows$status == "stale"))
    expect_true(closeProjectRegistry(inspection))

    for (context in contexts[artifact_names]) {
        for (path_name in c("staging", "trash")) {
            path <- artifactPath(context$getPaths(), path_name)
            entries <- if (dir.exists(path)) {
                list.files(path, all.files = TRUE, no.. = TRUE)
            } else {
                character()
            }
            expect_length(entries, 0L)
        }
    }
    expect_identical(dia_memory_manager$getState(), dia_object)
    expect_identical(nondia_memory_manager$getState(), nondia_object)

    globals_after <- vapply(
        globals,
        exists,
        logical(1),
        envir = .GlobalEnv,
        inherits = FALSE
    )
    expect_identical(globals_after, globals_before)
})

test_that("legacy artifact roots retain their public registry workflow ID", {
    for (package in c("DBI", "duckdb", "filelock", "arrow")) {
        skip_if_not_installed(package)
    }
    root <- withr::local_tempdir(pattern = "omics-art-047-legacy-")
    context <- omics014Context(root, "hybrid-legacy-project")
    identity <- context$getIdentity()
    store <- newArtifactStore(context$getPaths(), identity$project_id)
    expect_true(artifactWorkflowStateEnsureRootManifest(store, identity))
    expect_false(file.exists(artifactStoreResolveFile(
        store,
        artifactRegistryScopePath(store)
    )))

    object <- omics014Protein()
    manager <- omics014Manager(context)
    manager$saveState(
        "legacy_dia_state",
        object,
        object@args,
        "legacy registry compatibility"
    )
    generation_id <- manager$getCurrentGenerationId()
    expect_true(manager$close())
    hybridExpectClosed(manager)
    expect_identical(
        artifactRegistryIdentity(store, identity)$workflow_id,
        identity$workflow_id
    )

    reopened <- omics014Manager(context)
    expect_identical(reopened$getCurrentGenerationId(), generation_id)
    expect_identical(reopened$getState(), object)
    expect_true(reopened$close())
    hybridExpectClosed(reopened)

    inspection <- openProjectRegistryReadOnly(projectRegistryForContext(context))
    workflows <- projectRegistryQuery(inspection, "workflows")
    expect_identical(workflows$workflow_id, identity$workflow_id)
    expect_true(closeProjectRegistry(inspection))
})

test_that("registry scope metadata is immutable and fail closed", {
    for (package in c("DBI", "duckdb", "filelock", "arrow")) {
        skip_if_not_installed(package)
    }
    root <- withr::local_tempdir(pattern = "omics-art-047-scope-")
    context <- omics014Context(root, "hybrid-scope-project")
    manager <- omics014Manager(context)
    expect_true(manager$close())
    hybridExpectClosed(manager)

    identity <- context$getIdentity()
    store <- newArtifactStore(context$getPaths(), identity$project_id)
    scope_path <- artifactStoreResolveFile(
        store,
        artifactRegistryScopePath(store),
        must_exist = TRUE
    )
    scope <- jsonlite::read_json(scope_path, simplifyVector = TRUE)
    scope$registry_workflow_id <- paste0(
        scope$registry_workflow_id,
        "-tampered"
    )
    jsonlite::write_json(
        scope,
        scope_path,
        auto_unbox = TRUE,
        pretty = TRUE
    )
    expect_error(
        omics014Manager(context),
        class = "multischolar_incompatible_artifact_registry_scope"
    )
})

test_that("hybrid evidence remains public operational and tuple specific", {
    repo <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
    paths <- c(
        metabolomics = file.path(
            repo,
            "tests/testdata/omics-parity/metabolomics/workloads/",
            "mixed-public-representative-v1.json"
        ),
        lipidomics = file.path(
            repo,
            "tests/testdata/omics-parity/lipidomics/workloads/",
            "mixed-public-representative-v1.json"
        )
    )
    contracts <- lapply(paths, jsonlite::read_json, simplifyVector = FALSE)
    expect_true(all(vapply(
        contracts,
        \(contract) identical(contract$privacy$classification, "public_synthetic"),
        logical(1)
    )))
    expect_identical(anyDuplicated(vapply(
        contracts,
        `[[`,
        character(1),
        "workload_id"
    )), 0L)
    expect_false(any(grepl(
        "private|cotton|path|fasta|sequence|header",
        workflowSessionJson(contracts),
        ignore.case = TRUE
    )))
    expect_true(all(vapply(
        artifactMetabolomicsCloseoutDecisions(),
        \(decision) identical(decision$promotion_status, "withheld"),
        logical(1)
    )))
    expect_true(all(vapply(
        artifactLipidomicsCloseoutDecisions(),
        \(decision) identical(decision$promotion_status, "withheld"),
        logical(1)
    )))
})
