# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

.PROT_NONDIA_SESSION_IMPORT_ROLES <- "canonical_data"
.PROT_NONDIA_SESSION_DESIGN_ROLES <- c(
    "cleaned_data", "design_matrix", "contrasts", "args", "annotations",
    "sequences"
)

#' Resolve the independent non-DIA artifact session switch
#' @param direction Session export or restore direction.
#' @return One of `"enabled"` or `"disabled"`.
#' @noRd
protNonDiaSessionMode <- function(direction = c("export", "restore")) {
    direction <- match.arg(direction)
    option <- paste0(
        "multischolar.prot_nondia.artifact_session_",
        direction
    )
    match.arg(getOption(option, "enabled"), c("enabled", "disabled"))
}

#' Resolve non-DIA session manager resources
#' @param workflow_data Mutable proteomics workflow state.
#' @return Exact descriptor, manager, and ownership flag, or `NULL`.
#' @noRd
protNonDiaSessionManagerResources <- function(workflow_data) {
    descriptor <- protNonDiaNormDescriptor(workflow_data)
    if (is.null(descriptor)) return(NULL)
    manager <- workflow_data$state_manager
    if (inherits(manager, "ArtifactWorkflowState")) {
        return(list(
            descriptor = descriptor,
            manager = manager,
            owned = FALSE
        ))
    }
    artifact_manager <- protNonDiaNormArtifactManager(
        workflow_data,
        manager,
        descriptor
    )$manager
    valid <- identical(
        workflowStateCurrentName(manager),
        artifact_manager$getCurrentStateName()
    ) && identical(
        workflowStateHistory(manager),
        artifact_manager$getHistory()
    ) && identical(manager$getState(), artifact_manager$getState())
    if (!isTRUE(valid)) {
        artifact_manager$close()
        workflowSessionAbort(
            "non-DIA session managers differ at export",
            "multischolar_prot_nondia_session_manager_mismatch"
        )
    }
    list(
        descriptor = descriptor,
        manager = artifact_manager,
        owned = TRUE
    )
}

#' Test exact non-DIA artifact session eligibility
#' @param workflow_data Mutable proteomics workflow state.
#' @param direction Session export or restore direction.
#' @return A scalar logical.
#' @noRd
protNonDiaSessionArtifactEligible <- function(
    workflow_data,
    direction = "export"
) {
    if (!identical(protNonDiaSessionMode(direction), "enabled") ||
        is.null(workflow_data)) {
        return(FALSE)
    }
    descriptor <- protNonDiaNormDescriptor(workflow_data)
    if (is.null(descriptor)) return(FALSE)
    result <- workflow_data$artifact_stage_results$correlation_filter
    manager <- workflow_data$state_manager
    current <- tryCatch(workflowStateCurrentName(manager), error = \(...) NULL)
    artifact_ready <- inherits(manager, "ArtifactWorkflowState") ||
        (is.list(result) && isTRUE(result$committed))
    artifact_ready && identical(current, "correlation_filtered")
}

#' Select bounded portable non-DIA session metadata
#' @param session_data Collected normalization session data.
#' @return A runtime-free metadata list.
#' @noRd
protNonDiaSessionPortableMetadata <- function(session_data) {
    fields <- c(
        .PROT_DIA_SESSION_METADATA_FIELDS,
        "normalization_state_refs", "final_for_da_ref"
    )
    metadata <- session_data[intersect(fields, names(session_data))]
    metadata$ruv_optimization_result <- protDiaSessionPortableRuvResult(
        metadata$ruv_optimization_result
    )
    metadata
}

#' Resolve one non-DIA import or design stage result
#' @param workflow_data Mutable proteomics workflow state.
#' @param stage_id Import or design stage identifier.
#' @return Run, generation, and normalized refs.
#' @noRd
protNonDiaSessionStageResult <- function(workflow_data, stage_id) {
    result <- workflow_data$artifact_stage_results[[stage_id]]
    refs <- result$refs
    if (!is.list(refs) || length(refs) == 0L) {
        refs <- workflow_data$artifact_readthrough_refs[[stage_id]]
    }
    if (!is.list(refs) || length(refs) == 0L) {
        workflowSessionAbort(
            sprintf("non-DIA session is missing %s refs", stage_id),
            "multischolar_invalid_prot_nondia_session_stage"
        )
    }
    refs <- lapply(refs, artifactStoreNormalizeRef)
    run_ids <- unique(unlist(lapply(refs, `[[`, "provenance_ids")))
    generations <- unique(vapply(
        refs,
        \(ref) ref$logical_key$generation_id,
        character(1)
    ))
    run_id <- result$run_id
    generation_id <- result$generation_id
    proof <- workflow_data$artifact_readthrough_proof
    proof_run_id <- if (identical(stage_id, "import")) {
        proof$import_run_id
    } else {
        proof$design_run_id
    }
    if (!workflowCapabilityScalarString(run_id) &&
        workflowCapabilityScalarString(proof_run_id)) {
        run_id <- proof_run_id
    }
    if (!workflowCapabilityScalarString(run_id) && length(run_ids) == 1L) {
        run_id <- run_ids[[1L]]
    }
    if (!workflowCapabilityScalarString(generation_id) &&
        length(generations) == 1L) {
        generation_id <- generations[[1L]]
    }
    list(run_id = run_id, generation_id = generation_id, refs = refs)
}

#' Build parent-linked non-DIA stage evidence for a session
#' @param workflow_data Mutable proteomics workflow state.
#' @param descriptor Exact supported workflow descriptor.
#' @return Import and design session evidence.
#' @noRd
protNonDiaSessionStageEvidence <- function(workflow_data, descriptor) {
    import <- protNonDiaSessionStageResult(workflow_data, "import")
    design <- protNonDiaSessionStageResult(workflow_data, "design")
    list(
        import = list(
            run_id = import$run_id,
            generation_id = import$generation_id,
            parameters = list(
                capability_id = descriptor$descriptor_id,
                readthrough_contract_version = 1L,
                column_mapping_serialized =
                    artifactWorkflowStateSerializeMetadata(
                        workflow_data$column_mapping,
                        "non-DIA session column mapping"
                    )
            ),
            refs = import$refs
        ),
        design = list(
            run_id = design$run_id,
            generation_id = design$generation_id,
            parameters = list(
                capability_id = descriptor$descriptor_id,
                readthrough_contract_version = 1L,
                contrasts_kind = artifactStageContrastsKind(
                    workflow_data$contrasts_tbl
                ),
                annotation_available = !is.null(
                    workflow_data$uniprot_dat_cln
                ),
                sequence_data_available = !is.null(
                    workflow_data$aa_seq_tbl_final
                )
            ),
            refs = design$refs
        )
    )
}

#' Build a bounded non-DIA artifact session manifest
#' @param workflow_data Mutable proteomics workflow state.
#' @param norm_data Mutable normalization module state.
#' @param session_data Collected legacy-compatible session data.
#' @param compatibility_path Current compatibility RDS path.
#' @return A validated workflow session manifest, or `NULL` when ineligible.
#' @noRd
buildProtNonDiaSessionManifest <- function(
    workflow_data,
    norm_data,
    session_data,
    compatibility_path
) {
    if (!protNonDiaSessionArtifactEligible(workflow_data, "export")) {
        return(NULL)
    }
    resources <- protNonDiaSessionManagerResources(workflow_data)
    if (is.null(resources)) return(NULL)
    if (isTRUE(resources$owned)) on.exit(resources$manager$close(), add = TRUE)
    context <- workflow_data$workflow_context
    manager <- resources$manager
    descriptor <- resources$descriptor
    identity <- context$getIdentity()
    store <- newArtifactStore(context$getPaths(), identity$project_id)
    snapshot <- workflowStateManifest(manager)
    lineage <- protDiaSessionManifestForLineage(store, snapshot)
    current_manifest <- artifactWorkflowStateReadManifest(
        store,
        tail(lineage, 1L)[[1L]]$manifest_relative_path
    )
    project_root <- context$getProjectRoot()
    manifest <- list(
        schema = .WORKFLOW_SESSION_SCHEMA,
        schema_version = .WORKFLOW_SESSION_VERSION,
        backend = "artifact",
        identity = identity[c(
            "project_id", "workflow_id", "omic_type", "omic_label",
            "workflow_slug"
        )],
        descriptor = list(
            descriptor_id = descriptor$descriptor_id,
            descriptor_version = descriptor$descriptor_version,
            descriptor_digest = descriptor$descriptor_digest
        ),
        workflow_state = list(
            schema = snapshot$schema,
            schema_version = snapshot$schema_version,
            backend = snapshot$backend,
            current_generation_id = snapshot$current_generation_id,
            current_state = snapshot$current_state,
            active_lineage = lineage,
            workflow_type = snapshot$workflow_type,
            audit_enabled = snapshot$audit_enabled,
            state_semantic_digest = current_manifest$data$semantic_digest
        ),
        stage_refs = protNonDiaSessionStageEvidence(
            workflow_data,
            descriptor
        ),
        plot_refs = protDiaSessionPlotRefs(norm_data, project_root),
        metadata_json = workflowSessionEncodeMetadata(
            protNonDiaSessionPortableMetadata(session_data),
            "non-DIA session metadata"
        ),
        compatibility = list(
            relative_path = workflowSessionProjectRelativePath(
                compatibility_path,
                project_root
            ),
            byte_digest = artifactByteDigest(compatibility_path),
            generation_id = snapshot$current_generation_id,
            state_semantic_digest = current_manifest$data$semantic_digest
        ),
        exported_at = format(
            session_data$export_timestamp,
            "%Y-%m-%dT%H:%M:%OS6Z",
            tz = "UTC"
        ),
        manifest_digest = NULL
    )
    manifest$manifest_digest <- workflowSessionContentDigest(manifest)
    validateWorkflowSessionManifest(manifest)
}

#' Publish immutable and latest non-DIA session manifests
#' @param workflow_data Mutable proteomics workflow state.
#' @param norm_data Mutable normalization module state.
#' @param session_data Collected legacy-compatible session data.
#' @param export_artifacts Legacy export paths.
#' @param source_dir Project source directory.
#' @param failure_injector Optional session failure injector used by tests.
#' @return Published session manifest paths and data.
#' @noRd
saveProtNonDiaSessionManifest <- function(
    workflow_data,
    norm_data,
    session_data,
    export_artifacts,
    source_dir,
    failure_injector = NULL
) {
    if (!protNonDiaSessionArtifactEligible(workflow_data, "export")) {
        return(list(enabled = FALSE, reason = "artifact_session_not_eligible"))
    }
    timestamp <- sub(
        "^filtered_session_data_(.*)[.]rds$",
        "\\1",
        export_artifacts$sessionFilename
    )
    manifest <- buildProtNonDiaSessionManifest(
        workflow_data,
        norm_data,
        session_data,
        export_artifacts$sessionFilepath
    )
    timestamped <- file.path(
        source_dir,
        sprintf("filtered_session_artifact_%s.json", timestamp)
    )
    latest <- file.path(source_dir, "filtered_session_artifact_latest.json")
    workflowSessionPublishImmutable(manifest, timestamped, failure_injector)
    workflowSessionReplaceLatest(manifest, latest, failure_injector)
    list(
        enabled = TRUE,
        ok = TRUE,
        reason = "artifact_session_published",
        manifest = manifest,
        sessionFilepath = timestamped,
        latestFilepath = latest
    )
}

#' Validate non-DIA session identity and descriptor
#' @param manifest Validated session manifest.
#' @param context Exact bound workflow context.
#' @param descriptor Exact supported workflow descriptor.
#' @return `TRUE`, invisibly.
#' @noRd
protNonDiaSessionValidateIdentity <- function(manifest, context, descriptor) {
    identity <- context$getIdentity()
    expected <- identity[names(manifest$identity)]
    expected_descriptor <- list(
        descriptor_id = descriptor$descriptor_id,
        descriptor_version = descriptor$descriptor_version,
        descriptor_digest = descriptor$descriptor_digest
    )
    if (!identical(manifest$identity, expected) ||
        !identical(manifest$descriptor, expected_descriptor)) {
        workflowSessionAbort(
            "non-DIA session identity or descriptor differs from the project",
            "multischolar_prot_nondia_session_identity_mismatch"
        )
    }
    invisible(TRUE)
}

#' Validate non-DIA session import/design refs against the registry
#' @param manifest Validated session manifest.
#' @param context Exact bound workflow context.
#' @param descriptor Exact supported workflow descriptor.
#' @return `TRUE`, invisibly.
#' @noRd
protNonDiaSessionValidateRegistry <- function(manifest, context, descriptor) {
    evidence <- collectProtNonDiaResumeEvidence(
        context,
        descriptor,
        payload_validation = "digest"
    )
    for (stage_id in c("import", "design")) {
        expected <- manifest$stage_refs[[stage_id]]
        actual <- evidence[[stage_id]]
        expected_refs <- lapply(expected$refs, artifactStoreNormalizeRef)
        actual_refs <- lapply(actual$refs, artifactStoreNormalizeRef)
        valid <- identical(expected$run_id, actual$run_id) &&
            identical(expected$generation_id, actual$generation_id) &&
            identical(expected_refs, actual_refs)
        if (!isTRUE(valid)) {
            workflowSessionAbort(
                sprintf(
                    "non-DIA session %s refs differ from the registry",
                    stage_id
                ),
                "multischolar_prot_nondia_session_registry_mismatch"
            )
        }
    }
    invisible(TRUE)
}

#' Validate one non-DIA session stage and immutable refs
#' @param store Validated artifact store.
#' @param stage Session stage evidence.
#' @param roles Exact expected roles.
#' @param stage_id Stage identifier.
#' @param adapter Exact descriptor read-through adapter.
#' @return Normalized stage evidence.
#' @noRd
protNonDiaSessionValidateStage <- function(
    store,
    stage,
    roles,
    stage_id,
    adapter
) {
    required <- c("run_id", "generation_id", "parameters", "refs")
    valid <- is.list(stage) && identical(names(stage), required) &&
        workflowCapabilityScalarString(stage$run_id) &&
        workflowCapabilityScalarString(stage$generation_id) &&
        is.list(stage$parameters) && is.list(stage$refs) &&
        identical(names(stage$refs), roles)
    if (!isTRUE(valid)) {
        workflowSessionAbort(
            sprintf("non-DIA session %s evidence is incomplete", stage_id),
            "multischolar_invalid_prot_nondia_session_stage"
        )
    }
    refs <- lapply(stage$refs, \(ref) {
        validated <- artifactStageValidateStoredRef(adapter, store, ref)$ref
        lineage_valid <- identical(
            validated$logical_key$stage_id,
            stage_id
        ) && identical(
            validated$logical_key$generation_id,
            stage$generation_id
        ) && stage$run_id %in% validated$provenance_ids
        if (!isTRUE(lineage_valid)) {
            workflowSessionAbort(
                sprintf("non-DIA session %s ref lineage is invalid", stage_id),
                "multischolar_invalid_prot_nondia_session_stage"
            )
        }
        validated
    })
    names(refs) <- roles
    stage$refs <- refs
    stage
}

#' Restore non-DIA design, contrasts, annotations, and sequences
#' @param store Validated artifact store.
#' @param stage_refs Validated import/design session stages.
#' @param adapter Exact descriptor read-through adapter.
#' @return Restored scientific workflow values.
#' @noRd
protNonDiaSessionRestoreScientificValues <- function(
    store,
    stage_refs,
    adapter
) {
    design <- stage_refs$design
    tables <- readArtifactStageRoles(
        adapter,
        store,
        design,
        c("design_matrix", "contrasts", "annotations", "sequences")
    )
    parameters <- design$parameters
    list(
        design_matrix = tables$design_matrix,
        contrasts_tbl = artifactStageRestoreContrasts(
            adapter,
            tables$contrasts,
            parameters$contrasts_kind
        ),
        uniprot_dat_cln = artifactStageRestoreOptional(
            adapter,
            tables$annotations,
            parameters$annotation_available,
            "annotations"
        ),
        aa_seq_tbl_final = artifactStageRestoreOptional(
            adapter,
            tables$sequences,
            parameters$sequence_data_available,
            "sequences"
        )
    )
}

#' Validate a non-DIA session compatibility fingerprint
#' @param manifest Validated session manifest.
#' @param project_root Moved project root.
#' @return `TRUE`, invisibly.
#' @noRd
protNonDiaSessionValidateCompatibility <- function(manifest, project_root) {
    fingerprint <- manifest$compatibility
    path <- artifactResolveContainedPath(
        project_root,
        fingerprint$relative_path,
        must_exist = TRUE
    )
    valid <- identical(artifactByteDigest(path), fingerprint$byte_digest) &&
        identical(
            fingerprint$generation_id,
            manifest$workflow_state$current_generation_id
        ) && identical(
            fingerprint$state_semantic_digest,
            manifest$workflow_state$state_semantic_digest
        )
    if (!isTRUE(valid)) {
        workflowSessionAbort(
            "non-DIA compatibility session fingerprint differs",
            "multischolar_prot_nondia_session_compatibility_mismatch"
        )
    }
    invisible(TRUE)
}

#' Validate project-contained non-DIA session plot fingerprints
#' @param plot_refs Named relative plot fingerprints.
#' @param project_root Moved project root.
#' @return `TRUE`, invisibly.
#' @noRd
protNonDiaSessionValidatePlots <- function(plot_refs, project_root) {
    for (ref in plot_refs) {
        required <- c("relative_path", "byte_digest")
        if (!is.list(ref) || !identical(names(ref), required)) {
            workflowSessionAbort(
                "non-DIA session plot ref is invalid",
                "multischolar_invalid_prot_nondia_session_plot"
            )
        }
        path <- artifactResolveContainedPath(
            project_root,
            ref$relative_path,
            must_exist = TRUE
        )
        if (!identical(artifactByteDigest(path), ref$byte_digest)) {
            workflowSessionAbort(
                "non-DIA session plot digest does not match",
                "multischolar_prot_nondia_session_plot_mismatch"
            )
        }
    }
    invisible(TRUE)
}

#' Restore a moved non-DIA artifact session manifest
#' @param manifest_path Session manifest path.
#' @param experiment_paths Moved workflow project paths.
#' @param resource_policy Optional project registry resource policy.
#' @return A verified artifact session bundle with an owned manager.
#' @noRd
restoreProtNonDiaSessionManifest <- function(
    manifest_path,
    experiment_paths,
    resource_policy = NULL
) {
    manifest <- readWorkflowSessionManifest(manifest_path)
    descriptor <- protNonDiaReadthroughDescriptor(
        manifest$descriptor$descriptor_id
    )
    paths <- experiment_paths
    paths$project_id <- manifest$identity$project_id
    paths$omic_label <- manifest$identity$omic_label
    prepared <- createProtNonDiaResumeContext(
        paths,
        manifest$identity$omic_label,
        descriptor$descriptor_id,
        storage_policy = list(
            requested_backend = "artifact",
            requested_rollout = "dual_write",
            migration_requested = TRUE,
            project_id = manifest$identity$project_id
        )
    )
    if (!isTRUE(prepared$enabled)) {
        workflowSessionAbort(
            "non-DIA artifact session cannot bind the moved project",
            "multischolar_prot_nondia_session_context_unavailable"
        )
    }
    context <- prepared$context
    protNonDiaSessionValidateIdentity(manifest, context, descriptor)
    protNonDiaSessionValidateRegistry(manifest, context, descriptor)
    manager <- newProtNonDiaResumeStateManager(
        context,
        descriptor,
        resource_policy
    )
    manager_owned <- TRUE
    on.exit(if (manager_owned) try(manager$close(), silent = TRUE), add = TRUE)
    identity <- context$getIdentity()
    store <- newArtifactStore(context$getPaths(), identity$project_id)
    current_manifest <- protDiaSessionValidateLineage(
        manifest,
        manager,
        store
    )
    adapter <- protNonDiaArtifactReadthroughAdapter(descriptor)
    stages <- list(
        import = protNonDiaSessionValidateStage(
            store,
            manifest$stage_refs$import,
            .PROT_NONDIA_SESSION_IMPORT_ROLES,
            "import",
            adapter
        ),
        design = protNonDiaSessionValidateStage(
            store,
            manifest$stage_refs$design,
            .PROT_NONDIA_SESSION_DESIGN_ROLES,
            "design",
            adapter
        )
    )
    scientific <- protNonDiaSessionRestoreScientificValues(
        store,
        stages,
        adapter
    )
    state_object <- manager$getState()
    role <- artifactProteomicsNonDiaCodecRole(
        descriptor$descriptor_id,
        manager$getCurrentStateName(),
        descriptor$descriptor_version
    )
    artifactValidateProteomicsNonDiaProteinState(
        state_object,
        role,
        manager$getCurrentStateName()
    )
    valid <- identical(
        methods::validObject(state_object, test = TRUE),
        TRUE
    ) && identical(
        current_manifest$data$semantic_digest,
        manifest$workflow_state$state_semantic_digest
    ) && identical(state_object@design_matrix, scientific$design_matrix)
    if (!isTRUE(valid)) {
        workflowSessionAbort(
            "non-DIA session scientific state failed exact validation",
            "multischolar_prot_nondia_session_scientific_mismatch"
        )
    }
    protNonDiaSessionValidatePlots(
        manifest$plot_refs,
        context$getProjectRoot()
    )
    protNonDiaSessionValidateCompatibility(manifest, context$getProjectRoot())
    metadata <- workflowSessionDecodeMetadata(
        manifest$metadata_json,
        "non-DIA session metadata"
    )
    session_data <- c(list(
        r6_current_state_name = manager$getCurrentStateName(),
        r6_complete_states = manager$states,
        r6_state_history = manager$getHistory(),
        current_s4_object = state_object,
        correlation_filtered_s4 = state_object,
        contrasts_tbl = scientific$contrasts_tbl,
        design_matrix = scientific$design_matrix,
        config_list = state_object@args
    ), metadata)
    session_data$workflow_state_manifest <- workflowStateManifest(manager)
    session_data$artifact_session_manifest <- manifest
    session_data$uniprot_dat_cln <- scientific$uniprot_dat_cln
    session_data$aa_seq_tbl_final <- scientific$aa_seq_tbl_final
    validateProtDaFilteredSession(
        session_data,
        source_dir = experiment_paths$source_dir
    )
    manager_owned <- FALSE
    list(
        format = "artifact",
        session_data = session_data,
        state_manager = manager,
        workflow_context = context,
        uniprot_dat_cln = scientific$uniprot_dat_cln,
        aa_seq_tbl_final = scientific$aa_seq_tbl_final,
        manifest = manifest
    )
}

#' Reconstruct a current legacy session from a non-DIA artifact bundle
#' @param bundle Verified non-DIA artifact session bundle.
#' @return A legacy-compatible current session payload.
#' @noRd
reconstructProtNonDiaCompatibilitySession <- function(bundle) {
    if (!is.list(bundle) || !identical(bundle$format, "artifact") ||
        !inherits(bundle$state_manager, "ArtifactWorkflowState")) {
        workflowSessionAbort(
            "non-DIA compatibility reconstruction requires an artifact session",
            "multischolar_invalid_prot_nondia_session_reconstruction"
        )
    }
    snapshot <- materializeProtDiaLegacyStateSnapshot(bundle$state_manager)
    session_data <- bundle$session_data
    session_data$r6_current_state_name <- snapshot$r6_current_state_name
    session_data$r6_complete_states <- snapshot$r6_complete_states
    session_data$r6_state_history <- snapshot$r6_state_history
    session_data$workflow_state_manifest <- NULL
    session_data$artifact_session_manifest <- NULL
    session_data$artifact_source_generation <-
        bundle$manifest$workflow_state$current_generation_id
    session_data$artifact_source_state_digest <-
        bundle$manifest$workflow_state$state_semantic_digest
    session_data
}

#' Atomically write a current non-DIA legacy compatibility session
#' @param bundle Verified non-DIA artifact session bundle.
#' @param path Destination RDS path.
#' @param save_rds_fn Injectable RDS writer used by tests.
#' @return Compatibility path and integrity fingerprint.
#' @noRd
writeProtNonDiaCompatibilitySession <- function(
    bundle,
    path,
    save_rds_fn = saveRDS
) {
    session_data <- reconstructProtNonDiaCompatibilitySession(bundle)
    temporary <- file.path(
        dirname(path),
        paste0(".", basename(path), ".", artifactOpaqueId("tmp"), ".tmp")
    )
    on.exit(if (file.exists(temporary)) unlink(temporary, force = FALSE), add = TRUE)
    save_rds_fn(session_data, temporary)
    validateProtDaFilteredSession(readRDS(temporary), dirname(path))
    if (!isTRUE(file.rename(temporary, path))) {
        workflowSessionAbort(
            "non-DIA compatibility session could not be atomically published",
            "multischolar_prot_nondia_session_compatibility_write_failed"
        )
    }
    list(
        path = path,
        byte_digest = artifactByteDigest(path),
        generation_id = bundle$manifest$workflow_state$current_generation_id,
        state_semantic_digest =
            bundle$manifest$workflow_state$state_semantic_digest
    )
}
