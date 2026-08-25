# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

metabSessionMode <- function(direction = c("export", "restore")) {
    direction <- match.arg(direction)
    option <- paste0("multischolar.metabolomics.artifact_session_", direction)
    match.arg(getOption(option, "enabled"), c("enabled", "disabled"))
}

metabSessionManagerResources <- function(workflow_data) {
    manager <- workflow_data$state_manager
    if (inherits(manager, "ArtifactWorkflowState")) {
        return(list(manager = manager, owned = FALSE))
    }
    if (!metabQcWorkflowData(workflow_data) ||
        !metabArtifactCoordinatorOwned(workflow_data)) {
        return(NULL)
    }
    resources <- metabQcArtifactManager(workflow_data, manager)
    artifact_manager <- resources$manager
    valid <- identical(
        metabQcCurrentStateName(manager),
        artifact_manager$getCurrentStateName()
    ) && identical(manager$getState(), artifact_manager$getState())
    if (!isTRUE(valid)) {
        artifact_manager$close()
        workflowSessionAbort(
            "metabolomics session managers differ at export",
            "multischolar_metab_session_manager_mismatch"
        )
    }
    list(manager = artifact_manager, owned = TRUE)
}

metabSessionArtifactEligible <- function(workflow_data, direction = "export") {
    if (!identical(metabSessionMode(direction), "enabled") ||
        !metabQcWorkflowData(workflow_data)) {
        return(FALSE)
    }
    manager <- workflow_data$state_manager
    current <- tryCatch(metabQcCurrentStateName(manager), error = function(...) NULL)
    ready_state <- current %in% c(
        "metab_norm_complete", "metab_correlation_filtered",
        "metab_ruv_corrected", "metab_normalized"
    )
    artifact_ready <- inherits(manager, "ArtifactWorkflowState") ||
        isTRUE(workflow_data$artifact_stage_results$correlation_skip$committed) ||
        isTRUE(workflow_data$artifact_stage_results$correlation_filter$committed)
    isTRUE(ready_state) && isTRUE(artifact_ready)
}

metabSessionPortableMetadata <- function(session_data) {
    fields <- c(
        "normalization_method", "ruv_mode", "itsd_applied",
        "itsd_aggregation", "log_offset", "correlation_threshold",
        "ruv_grouping_variable", "feature_counts", "metabolite_counts",
        "assay_names", "experiment_label", "omic_type", "contrasts_tbl"
    )
    session_data[intersect(fields, names(session_data))]
}

metabSessionLineage <- function(store, snapshot) {
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

metabSessionStateRefs <- function(manager) {
    rows <- manager$states
    lapply(names(rows), function(state_name) {
        row <- rows[[state_name]]
        list(
            logical_name = state_name,
            generation_id = row$generation_id,
            manifest_relative_path = row$manifest_relative_path
        )
    })
}

buildMetabSessionManifest <- function(
    workflow_data,
    session_data,
    compatibility_path
) {
    if (!metabSessionArtifactEligible(workflow_data, "export")) return(NULL)
    resources <- metabSessionManagerResources(workflow_data)
    if (is.null(resources)) return(NULL)
    if (isTRUE(resources$owned)) on.exit(resources$manager$close(), add = TRUE)
    manager <- resources$manager
    context <- workflow_data$workflow_context
    descriptor <- artifactMetabolomicsWorkflowDescriptor()
    identity <- context$getIdentity()
    store <- newArtifactStore(context$getPaths(), identity$project_id)
    snapshot <- workflowStateManifest(manager)
    lineage <- metabSessionLineage(store, snapshot)
    current_manifest <- artifactWorkflowStateReadManifest(
        store,
        tail(lineage, 1L)[[1L]]$manifest_relative_path
    )
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
        stage_refs = list(states = metabSessionStateRefs(manager)),
        plot_refs = list(),
        metadata_json = workflowSessionEncodeMetadata(
            metabSessionPortableMetadata(session_data),
            "metabolomics session metadata"
        ),
        compatibility = list(
            relative_path = workflowSessionProjectRelativePath(
                compatibility_path,
                context$getProjectRoot()
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

saveMetabSessionManifest <- function(
    workflow_data,
    session_data,
    export_files,
    source_dir,
    failure_injector = NULL
) {
    if (!metabSessionArtifactEligible(workflow_data, "export")) {
        return(list(enabled = FALSE, reason = "artifact_session_not_eligible"))
    }
    timestamp <- sub(
        "^metab_filtered_session_data_(.*)[.]rds$",
        "\\1",
        export_files$sessionFilename
    )
    manifest <- buildMetabSessionManifest(
        workflow_data,
        session_data,
        export_files$sessionFilepath
    )
    timestamped <- file.path(
        source_dir,
        sprintf("metab_filtered_session_artifact_%s.json", timestamp)
    )
    latest <- file.path(
        source_dir,
        "metab_filtered_session_artifact_latest.json"
    )
    workflowSessionPublishImmutable(manifest, timestamped, failure_injector)
    workflowSessionReplaceLatest(manifest, latest, failure_injector)
    result <- list(
        enabled = TRUE,
        ok = TRUE,
        manifest = manifest,
        sessionFilepath = timestamped,
        latestFilepath = latest
    )
    recordArtifactStageResult(workflow_data, "session_export", result)
    result
}

saveMetabSessionManifestSafely <- function(...) {
    args <- list(...)
    workflow_data <- args$workflow_data
    tryCatch(
        do.call(saveMetabSessionManifest, args),
        error = function(error) {
            result <- list(
                enabled = TRUE,
                ok = FALSE,
                reason = "artifact_session_export_failed",
                error_class = class(error),
                error_message = conditionMessage(error)
            )
            if (metabQcWorkflowData(workflow_data)) {
                recordArtifactStageResult(
                    workflow_data,
                    "session_export",
                    result
                )
            }
            result
        }
    )
}

metabSessionValidateIdentity <- function(manifest, context) {
    identity <- context$getIdentity()
    expected <- identity[names(manifest$identity)]
    descriptor <- artifactMetabolomicsWorkflowDescriptor()
    valid <- identical(manifest$identity, expected) &&
        identical(manifest$descriptor$descriptor_id, descriptor$descriptor_id) &&
        identical(
            manifest$descriptor$descriptor_version,
            descriptor$descriptor_version
        ) && identical(
            manifest$descriptor$descriptor_digest,
            descriptor$descriptor_digest
        )
    if (!isTRUE(valid)) {
        workflowSessionAbort(
            "metabolomics session identity or descriptor differs",
            "multischolar_metab_session_identity_mismatch"
        )
    }
    invisible(TRUE)
}

metabSessionValidateLineage <- function(manifest, manager, store) {
    current <- manager$getCurrentGenerationId()
    if (!identical(current, manifest$workflow_state$current_generation_id) ||
        !identical(
            manager$getCurrentStateName(),
            manifest$workflow_state$current_state
        )) {
        workflowSessionAbort(
            "metabolomics session current generation differs from project",
            "multischolar_metab_session_lineage_mismatch"
        )
    }
    for (entry in manifest$workflow_state$active_lineage) {
        state_manifest <- artifactWorkflowStateReadManifest(
            store,
            entry$manifest_relative_path
        )
        if (!identical(state_manifest$manifest_digest, entry$manifest_digest)) {
            workflowSessionAbort(
                "metabolomics session lineage fingerprint differs",
                "multischolar_metab_session_lineage_mismatch"
            )
        }
    }
    invisible(TRUE)
}

metabSessionCompatibilityStatus <- function(manifest, project_root) {
    path <- artifactResolveContainedPath(
        project_root,
        manifest$compatibility$relative_path
    )
    if (!file.exists(path)) {
        return(list(available = FALSE, valid = FALSE, path = path))
    }
    valid <- identical(
        artifactByteDigest(path),
        manifest$compatibility$byte_digest
    )
    if (!isTRUE(valid)) {
        workflowSessionAbort(
            "metabolomics compatibility RDS fingerprint differs",
            "multischolar_metab_session_compatibility_mismatch"
        )
    }
    list(available = TRUE, valid = TRUE, path = path)
}

buildMetabRestoredSessionData <- function(manifest, manager) {
    metadata <- workflowSessionDecodeMetadata(
        manifest$metadata_json,
        "metabolomics session metadata"
    )
    object <- manager$getState()
    c(metadata, list(
        r6_current_state_name = manager$getCurrentStateName(),
        current_s4_object = object,
        design_matrix = object@design_matrix,
        config_list = manager$getStateConfig(),
        assay_names = names(object@metabolite_data),
        normalization_complete = TRUE,
        ruv_complete = !is.null(object@args$ruvIII_C_Varying),
        correlation_filtering_complete = manager$getCurrentStateName() %in%
            c("metab_correlation_filtered", "metab_norm_complete")
    ))
}

restoreMetabSessionManifest <- function(path, context, manager = NULL) {
    if (!identical(metabSessionMode("restore"), "enabled")) {
        workflowSessionAbort(
            "metabolomics artifact session restore is disabled",
            "multischolar_metab_session_restore_disabled"
        )
    }
    manifest <- readWorkflowSessionManifest(path)
    metabSessionValidateIdentity(manifest, context)
    owned <- is.null(manager)
    if (owned) {
        manager <- ArtifactWorkflowState$new(
            workflow_context = context,
            dehydrate_fn = dehydrateMetabolomicsS4Artifact,
            validate_bundle_fn = validateMetabolomicsS4Bundle,
            hydrate_fn = hydrateMetabolomicsS4Artifact,
            descriptor_contract = artifactStageDescriptorContract(
                artifactMetabolomicsWorkflowDescriptor()
            )
        )
    }
    on.exit(if (owned) manager$close(), add = TRUE)
    store <- newArtifactStore(
        context$getPaths(),
        context$getIdentity()$project_id
    )
    metabSessionValidateLineage(manifest, manager, store)
    compatibility <- metabSessionCompatibilityStatus(
        manifest,
        context$getProjectRoot()
    )
    session_data <- buildMetabRestoredSessionData(manifest, manager)
    owned <- FALSE
    list(
        manifest = manifest,
        manager = manager,
        session_data = session_data,
        compatibility = compatibility
    )
}

writeMetabCompatibilitySession <- function(
    bundle,
    path,
    save_rds = saveRDS
) {
    object <- bundle$session_data$current_s4_object
    temporary <- paste0(path, ".", artifactOpaqueId("tmp"), ".tmp")
    on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
    save_rds(bundle$session_data, temporary)
    check <- readRDS(temporary)
    if (!identical(check$current_s4_object, object) ||
        !isTRUE(file.rename(temporary, path))) {
        workflowSessionAbort(
            "metabolomics compatibility RDS could not be verified",
            "multischolar_metab_session_reconstruction_failed"
        )
    }
    invisible(path)
}

readMetabArtifactOrLegacySession <- function(
    legacy_path,
    workflow_data,
    read_rds = readRDS
) {
    manifest_path <- file.path(
        dirname(legacy_path),
        "metab_filtered_session_artifact_latest.json"
    )
    context <- workflow_data$workflow_context
    if (identical(metabSessionMode("restore"), "enabled") &&
        file.exists(manifest_path) && inherits(context, "WorkflowContext") &&
        context$isBound() &&
        identical(context$getStorageDecision()$effective_backend, "artifact")) {
        restored <- tryCatch(
            restoreMetabSessionManifest(
                manifest_path,
                context,
                manager = if (inherits(
                    workflow_data$state_manager,
                    "ArtifactWorkflowState"
                )) {
                    workflow_data$state_manager
                } else {
                    NULL
                }
            ),
            error = function(error) error
        )
        if (!inherits(restored, "error")) {
            session_data <- restored$session_data
            if (!identical(restored$manager, workflow_data$state_manager)) {
                restored$manager$close()
            }
            return(session_data)
        }
    }
    if (!file.exists(legacy_path)) {
        workflowSessionAbort(
            "no valid metabolomics artifact or legacy session is available",
            "multischolar_missing_metab_session"
        )
    }
    read_rds(legacy_path)
}
