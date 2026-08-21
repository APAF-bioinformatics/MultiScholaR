# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

.PROT_DIA_SESSION_IMPORT_ROLES <- "canonical_data"
.PROT_DIA_SESSION_DESIGN_ROLES <- c(
    "cleaned_data", "design_matrix", "contrasts", "args", "annotations",
    "sequences"
)
.PROT_DIA_SESSION_METADATA_FIELDS <- c(
    "export_timestamp", "normalization_method", "ruv_mode", "ruv_applied",
    "ruv_k", "correlation_threshold", "workflow_type", "report_template",
    "use_limpa", "limpa_applied", "limpa_parameters",
    "limpa_dpc_quant_results", "fasta_metadata",
    "accession_cleanup_results", "ruv_optimization_result", "qc_params",
    "protein_counts", "mixed_species_analysis", "final_protein_count",
    "final_sample_count"
)

protDiaSessionMode <- function(direction = c("export", "restore")) {
    direction <- match.arg(direction)
    option <- paste0("multischolar.prot_dia.artifact_session_", direction)
    match.arg(getOption(option, "enabled"), c("enabled", "disabled"))
}

protDiaSessionArtifactEligible <- function(workflow_data, direction = "export") {
    if (!identical(protDiaSessionMode(direction), "enabled") ||
        is.null(workflow_data)) {
        return(FALSE)
    }
    manager <- tryCatch(
        workflow_data$state_manager,
        error = function(error) NULL
    )
    context <- tryCatch(
        workflow_data$workflow_context,
        error = function(error) NULL
    )
    if (!inherits(manager, "ArtifactWorkflowState") ||
        !inherits(context, "WorkflowContext") || !context$isBound()) {
        return(FALSE)
    }
    identity <- context$getIdentity()
    decision <- context$getStorageDecision()
    identical(decision$effective_backend, "artifact") &&
        identical(workflowStateType(manager), "DIA") &&
        protDiaArtifactIdentityMatches(identity)
}

materializeProtDiaLegacyStateSnapshot <- function(state_manager) {
    if (!inherits(state_manager, "ArtifactWorkflowState")) {
        return(workflowStateLegacySnapshot(state_manager))
    }
    history <- workflowStateHistory(state_manager)
    states <- lapply(history, function(state_name) {
        state_manager$getState(state_name)
    })
    names(states) <- history
    list(
        r6_current_state_name = workflowStateCurrentName(state_manager),
        r6_complete_states = states,
        r6_state_history = history
    )
}

protDiaSessionPortableRuvResult <- function(value) {
    if (!is.list(value)) return(value)
    value$best_cancor_plot <- NULL
    value
}

protDiaSessionPortableMetadata <- function(session_data) {
    metadata <- session_data[intersect(
        .PROT_DIA_SESSION_METADATA_FIELDS,
        names(session_data)
    )]
    metadata$ruv_optimization_result <- protDiaSessionPortableRuvResult(
        metadata$ruv_optimization_result
    )
    metadata
}

protDiaSessionPlotRefs <- function(norm_data, project_root) {
    paths <- unlist(norm_data$qc_plot_paths, use.names = TRUE)
    if (length(paths) == 0L) return(list())
    refs <- list()
    for (index in seq_along(paths)) {
        path <- paths[[index]]
        if (!workflowCapabilityScalarString(path) || !file.exists(path) ||
            dir.exists(path)) {
            next
        }
        role <- names(paths)[[index]]
        if (!workflowCapabilityScalarString(role)) {
            role <- sprintf("plot_%04d", index)
        }
        refs[[role]] <- list(
            relative_path = workflowSessionProjectRelativePath(
                path,
                project_root
            ),
            byte_digest = artifactByteDigest(path)
        )
    }
    refs
}

protDiaSessionManifestForLineage <- function(store, snapshot) {
    lapply(snapshot$active_lineage, function(entry) {
        manifest <- artifactWorkflowStateReadManifest(
            store,
            entry$manifest_relative_path
        )
        list(
            generation_id = entry$generation_id,
            logical_name = entry$logical_name,
            manifest_relative_path = entry$manifest_relative_path,
            manifest_digest = manifest$manifest_digest
        )
    })
}

protDiaSessionStageResult <- function(workflow_data, stage_id) {
    result <- workflow_data$artifact_stage_results[[stage_id]]
    refs <- result$refs
    if (!is.list(refs) || length(refs) == 0L) {
        refs <- workflow_data$artifact_readthrough_refs[[stage_id]]
    }
    if (!is.list(refs) || length(refs) == 0L) {
        workflowSessionAbort(
            sprintf("DIA-NN session is missing %s refs", stage_id),
            "multischolar_invalid_prot_dia_session_stage"
        )
    }
    refs <- lapply(refs, artifactStoreNormalizeRef)
    run_ids <- unique(unlist(lapply(refs, `[[`, "provenance_ids")))
    generations <- unique(vapply(
        refs,
        function(ref) ref$logical_key$generation_id,
        character(1)
    ))
    run_id <- result$run_id
    generation_id <- result$generation_id
    if (!workflowCapabilityScalarString(run_id) && length(run_ids) == 1L) {
        run_id <- run_ids[[1L]]
    }
    if (!workflowCapabilityScalarString(generation_id) &&
        length(generations) == 1L) {
        generation_id <- generations[[1L]]
    }
    list(run_id = run_id, generation_id = generation_id, refs = refs)
}

protDiaSessionStageEvidence <- function(workflow_data) {
    import <- protDiaSessionStageResult(workflow_data, "import")
    design <- protDiaSessionStageResult(workflow_data, "design")
    list(
        import = list(
            run_id = import$run_id,
            generation_id = import$generation_id,
            parameters = list(
                readthrough_contract_version = 1L,
                column_mapping_serialized = artifactWorkflowStateSerializeMetadata(
                    workflow_data$column_mapping,
                    "DIA-NN session column mapping"
                )
            ),
            refs = import$refs
        ),
        design = list(
            run_id = design$run_id,
            generation_id = design$generation_id,
            parameters = list(
                readthrough_contract_version = 1L,
                contrasts_kind = protDiaArtifactContrastsKind(
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

buildProtDiaSessionManifest <- function(
    workflow_data,
    norm_data,
    session_data,
    compatibility_path
) {
    if (!protDiaSessionArtifactEligible(workflow_data, "export")) {
        return(NULL)
    }
    context <- workflow_data$workflow_context
    manager <- workflow_data$state_manager
    identity <- context$getIdentity()
    descriptor <- artifactDiaWorkflowDescriptor()
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
        stage_refs = protDiaSessionStageEvidence(workflow_data),
        plot_refs = protDiaSessionPlotRefs(norm_data, project_root),
        metadata_json = workflowSessionEncodeMetadata(
            protDiaSessionPortableMetadata(session_data)
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

saveProtDiaSessionManifest <- function(
    workflow_data,
    norm_data,
    session_data,
    export_artifacts,
    source_dir,
    failure_injector = NULL
) {
    if (!protDiaSessionArtifactEligible(workflow_data, "export")) {
        return(list(enabled = FALSE, reason = "artifact_session_not_eligible"))
    }
    timestamp <- sub(
        "^filtered_session_data_(.*)[.]rds$",
        "\\1",
        export_artifacts$sessionFilename
    )
    manifest <- buildProtDiaSessionManifest(
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
    workflowSessionPublishImmutable(
        manifest,
        timestamped,
        failure_injector
    )
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

protDiaSessionValidateIdentity <- function(manifest, context) {
    identity <- context$getIdentity()
    expected <- identity[names(manifest$identity)]
    if (!identical(manifest$identity, expected)) {
        workflowSessionAbort(
            "DIA-NN session identity differs from the moved project",
            "multischolar_prot_dia_session_identity_mismatch"
        )
    }
    descriptor <- artifactDiaWorkflowDescriptor()
    expected_descriptor <- list(
        descriptor_id = descriptor$descriptor_id,
        descriptor_version = descriptor$descriptor_version,
        descriptor_digest = descriptor$descriptor_digest
    )
    if (!identical(manifest$descriptor, expected_descriptor)) {
        workflowSessionAbort(
            "DIA-NN session descriptor differs from the installed contract",
            "multischolar_prot_dia_session_descriptor_mismatch"
        )
    }
    invisible(TRUE)
}

protDiaSessionValidateRegistry <- function(manifest, context) {
    evidence <- collectProtDiaResumeEvidence(
        context,
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
                    "DIA-NN session %s refs differ from the registry",
                    stage_id
                ),
                "multischolar_prot_dia_session_registry_mismatch"
            )
        }
    }
    invisible(TRUE)
}

protDiaSessionValidateLineage <- function(manifest, manager, store) {
    snapshot <- workflowStateManifest(manager)
    expected <- manifest$workflow_state
    fields <- c(
        "schema", "schema_version", "backend", "current_generation_id",
        "current_state", "workflow_type", "audit_enabled"
    )
    if (!identical(snapshot[fields], expected[fields])) {
        workflowSessionAbort(
            "DIA-NN session does not select the current artifact generation",
            "multischolar_prot_dia_session_generation_mismatch"
        )
    }
    actual <- protDiaSessionManifestForLineage(store, snapshot)
    if (!identical(actual, expected$active_lineage)) {
        workflowSessionAbort(
            "DIA-NN session lineage differs from committed state",
            "multischolar_prot_dia_session_lineage_mismatch"
        )
    }
    current <- artifactWorkflowStateReadManifest(
        store,
        tail(actual, 1L)[[1L]]$manifest_relative_path
    )
    if (!identical(
        current$data$semantic_digest,
        expected$state_semantic_digest
    )) {
        workflowSessionAbort(
            "DIA-NN session state digest differs from committed state",
            "multischolar_prot_dia_session_state_mismatch"
        )
    }
    current
}

protDiaSessionValidateStage <- function(store, stage, roles, stage_id) {
    required <- c("run_id", "generation_id", "parameters", "refs")
    if (!is.list(stage) || !identical(names(stage), required) ||
        !workflowCapabilityScalarString(stage$run_id) ||
        !workflowCapabilityScalarString(stage$generation_id) ||
        !is.list(stage$parameters) || !is.list(stage$refs) ||
        !identical(names(stage$refs), roles)) {
        workflowSessionAbort(
            sprintf("DIA-NN session %s evidence is incomplete", stage_id),
            "multischolar_invalid_prot_dia_session_stage"
        )
    }
    refs <- lapply(stage$refs, function(ref) {
        validated <- protDiaArtifactValidateStoredRef(store, ref)
        normalized <- validated$ref
        if (!identical(normalized$logical_key$stage_id, stage_id) ||
            !identical(
                normalized$logical_key$generation_id,
                stage$generation_id
            ) || !stage$run_id %in% normalized$provenance_ids) {
            workflowSessionAbort(
                sprintf("DIA-NN session %s ref lineage is invalid", stage_id),
                "multischolar_invalid_prot_dia_session_stage"
            )
        }
        normalized
    })
    names(refs) <- roles
    stage$refs <- refs
    stage
}

protDiaSessionRestoreScientificValues <- function(store, stage_refs) {
    design <- stage_refs$design
    tables <- readProtDiaArtifactStageRoles(
        store,
        design,
        c("design_matrix", "contrasts", "annotations", "sequences")
    )
    parameters <- design$parameters
    list(
        design_matrix = tables$design_matrix,
        contrasts_tbl = protDiaArtifactRestoreContrasts(
            tables$contrasts,
            parameters$contrasts_kind
        ),
        uniprot_dat_cln = protDiaArtifactRestoreOptional(
            tables$annotations,
            parameters$annotation_available,
            "annotations"
        ),
        aa_seq_tbl_final = protDiaArtifactRestoreOptional(
            tables$sequences,
            parameters$sequence_data_available,
            "sequences"
        )
    )
}

protDiaSessionValidatePlots <- function(plot_refs, project_root) {
    for (ref in plot_refs) {
        required <- c("relative_path", "byte_digest")
        if (!is.list(ref) || !identical(names(ref), required)) {
            workflowSessionAbort(
                "DIA-NN session plot ref is invalid",
                "multischolar_invalid_prot_dia_session_plot"
            )
        }
        path <- artifactResolveContainedPath(
            project_root,
            ref$relative_path,
            must_exist = TRUE
        )
        if (!identical(artifactByteDigest(path), ref$byte_digest)) {
            workflowSessionAbort(
                "DIA-NN session plot digest does not match",
                "multischolar_prot_dia_session_plot_mismatch"
            )
        }
    }
    invisible(TRUE)
}

buildProtDiaRestoredSessionData <- function(
    manifest,
    manager,
    state_object,
    scientific_values
) {
    metadata <- workflowSessionDecodeMetadata(manifest$metadata_json)
    session_data <- c(list(
        r6_current_state_name = manager$getCurrentStateName(),
        r6_complete_states = manager$states,
        r6_state_history = manager$getHistory(),
        current_s4_object = state_object,
        correlation_filtered_s4 = state_object,
        contrasts_tbl = scientific_values$contrasts_tbl,
        design_matrix = scientific_values$design_matrix,
        config_list = manager$getStateConfig()
    ), metadata)
    session_data$workflow_state_manifest <- workflowStateManifest(manager)
    session_data$artifact_session_manifest <- manifest
    session_data$uniprot_dat_cln <- scientific_values$uniprot_dat_cln
    session_data$aa_seq_tbl_final <- scientific_values$aa_seq_tbl_final
    session_data
}

restoreProtDiaSessionManifest <- function(
    manifest_path,
    experiment_paths,
    resource_policy = NULL
) {
    manifest <- readWorkflowSessionManifest(manifest_path)
    paths <- experiment_paths
    paths$project_id <- manifest$identity$project_id
    paths$omic_label <- manifest$identity$omic_label
    prepared <- createProtDiaResumeContext(
        paths,
        manifest$identity$omic_label,
        storage_policy = list(
            requested_backend = "artifact",
            requested_rollout = "dual_write",
            migration_requested = TRUE,
            project_id = manifest$identity$project_id
        )
    )
    if (!isTRUE(prepared$enabled)) {
        workflowSessionAbort(
            "DIA-NN artifact session cannot bind the moved project",
            "multischolar_prot_dia_session_context_unavailable"
        )
    }
    context <- prepared$context
    protDiaSessionValidateIdentity(manifest, context)
    protDiaSessionValidateRegistry(manifest, context)
    manager <- newProtDiaArtifactStateManager(context, resource_policy)
    manager_owned <- TRUE
    on.exit(if (manager_owned) try(manager$close(), silent = TRUE), add = TRUE)
    identity <- context$getIdentity()
    store <- newArtifactStore(context$getPaths(), identity$project_id)
    current_manifest <- protDiaSessionValidateLineage(
        manifest,
        manager,
        store
    )
    stages <- list(
        import = protDiaSessionValidateStage(
            store,
            manifest$stage_refs$import,
            .PROT_DIA_SESSION_IMPORT_ROLES,
            "import"
        ),
        design = protDiaSessionValidateStage(
            store,
            manifest$stage_refs$design,
            .PROT_DIA_SESSION_DESIGN_ROLES,
            "design"
        )
    )
    scientific_values <- protDiaSessionRestoreScientificValues(store, stages)
    state_object <- manager$getState()
    valid <- methods::is(state_object, "ProteinQuantitativeData") &&
        identical(methods::validObject(state_object, test = TRUE), TRUE) &&
        identical(
            current_manifest$data$semantic_digest,
            manifest$workflow_state$state_semantic_digest
        ) && identical(state_object@design_matrix, scientific_values$design_matrix)
    if (!isTRUE(valid)) {
        workflowSessionAbort(
            "DIA-NN session scientific state failed exact validation",
            "multischolar_prot_dia_session_scientific_mismatch"
        )
    }
    protDiaSessionValidatePlots(manifest$plot_refs, context$getProjectRoot())
    session_data <- buildProtDiaRestoredSessionData(
        manifest,
        manager,
        state_object,
        scientific_values
    )
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
        uniprot_dat_cln = scientific_values$uniprot_dat_cln,
        aa_seq_tbl_final = scientific_values$aa_seq_tbl_final,
        manifest = manifest
    )
}

reconstructProtDiaCompatibilitySession <- function(bundle) {
    if (!is.list(bundle) || !identical(bundle$format, "artifact") ||
        !inherits(bundle$state_manager, "ArtifactWorkflowState")) {
        workflowSessionAbort(
            "compatibility reconstruction requires a verified artifact session",
            "multischolar_invalid_prot_dia_session_reconstruction"
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

writeProtDiaCompatibilitySession <- function(
    bundle,
    path,
    save_rds_fn = saveRDS
) {
    session_data <- reconstructProtDiaCompatibilitySession(bundle)
    temporary <- file.path(
        dirname(path),
        paste0(".", basename(path), ".", artifactOpaqueId("tmp"), ".tmp")
    )
    on.exit(if (file.exists(temporary)) unlink(temporary, force = FALSE), add = TRUE)
    save_rds_fn(session_data, temporary)
    validateProtDaFilteredSession(readRDS(temporary), dirname(path))
    if (!isTRUE(file.rename(temporary, path))) {
        workflowSessionAbort(
            "compatibility session could not be atomically published",
            "multischolar_prot_dia_session_compatibility_write_failed"
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
