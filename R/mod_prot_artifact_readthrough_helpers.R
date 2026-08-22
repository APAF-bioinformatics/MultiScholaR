# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

protDiaArtifactValidateStoredRef <- function(store, ref) {
    ref <- artifactStoreNormalizeRef(ref)
    sidecar_path <- artifactStoreManagedPaths(
        store,
        ref$logical_key,
        ref$artifact_id
    )$sidecar
    sidecar <- artifactStoreReadSidecar(
        store,
        sidecar_path,
        validate_payload = FALSE
    )
    if (!identical(artifactStoreNormalizeRef(sidecar$artifact_ref), ref)) {
        protDiaArtifactAbort(
            "DIA-NN artifact ref differs from its immutable sidecar",
            "multischolar_prot_dia_registry_ref_mismatch",
            artifact_id = ref$artifact_id
        )
    }
    list(ref = ref, sidecar = sidecar)
}

protDiaArtifactReadTable <- function(store, ref) {
    validated <- protDiaArtifactValidateStoredRef(store, ref)
    payload_path <- artifactStoreResolveFile(
        store,
        validated$ref$relative_path,
        must_exist = TRUE
    )
    payload <- arrow::read_parquet(
        payload_path,
        as_data_frame = FALSE
    )
    value <- decodeArtifactRectangular(
        payload,
        validated$sidecar$codec_metadata
    )
    validateArtifactRefPayload(
        validated$ref,
        store$project_root,
        list(
            kind = validated$sidecar$codec_metadata$kind,
            rows = nrow(value),
            columns = ncol(value),
            payloads = 1L,
            bytes = unname(as.numeric(file.info(payload_path)$size))
        )
    )
    value
}

readProtDiaArtifactStageRoles <- function(store, stage, roles) {
    missing <- setdiff(roles, names(stage$refs))
    if (length(missing) > 0L) {
        protDiaArtifactAbort(
            sprintf(
                "DIA-NN stage '%s' is missing role(s): %s",
                stage$stage_id,
                paste(missing, collapse = ", ")
            ),
            "multischolar_incomplete_prot_dia_readthrough_contract",
            stage_id = stage$stage_id,
            missing_roles = missing
        )
    }
    values <- lapply(
        stage$refs[roles],
        \(ref) protDiaArtifactReadTable(store, ref)
    )
    names(values) <- roles
    values
}

readProtDiaArtifactStage <- function(store, stage) {
    values <- lapply(stage$refs, \(ref) protDiaArtifactReadTable(store, ref))
    names(values) <- names(stage$refs)
    generations <- unique(vapply(
        stage$refs,
        \(ref) ref$logical_key$generation_id,
        character(1)
    ))
    if (length(generations) != 1L ||
        !identical(generations[[1L]], stage$generation_id)) {
        protDiaArtifactAbort(
            "DIA-NN stage hydration combined different generations",
            "multischolar_mixed_prot_dia_artifact_generation",
            stage_id = stage$stage_id,
            run_id = stage$run_id
        )
    }
    values
}

protDiaArtifactContractVersion <- function(value) {
    identical(
        workflowStateVersionValue(value),
        .PROT_DIA_READTHROUGH_VERSION
    )
}

protDiaArtifactValidateResumeContract <- function(evidence) {
    design <- evidence$design
    import <- evidence$import
    design_ref <- design$refs$cleaned_data
    import_ref <- import$refs$canonical_data
    parameters <- design$parameters
    valid <- protDiaArtifactContractVersion(
        parameters$readthrough_contract_version
    ) && protDiaArtifactContractVersion(
        import$parameters$readthrough_contract_version
    ) && identical(parameters$workflow_type, "DIA") &&
        workflowCapabilityScalarString(parameters$state_name) &&
        identical(parameters$parent_import_run_id, import$run_id) &&
        identical(parameters$parent_import_generation_id, import$generation_id) &&
        identical(parameters$parent_import_artifact_id, import_ref$artifact_id) &&
        identical(
            parameters$parent_import_semantic_digest,
            import_ref$hash_policy$semantic$digest
        ) && identical(
            evidence$current_state$logical_name,
            parameters$state_name
        ) && identical(
            design_ref$logical_key$generation_id,
            design$generation_id
        )
    if (!isTRUE(valid)) {
        protDiaArtifactAbort(
            "DIA-NN import, design, and state evidence is not one lineage",
            "multischolar_incomplete_prot_dia_readthrough_contract",
            project_id = evidence$identity$project_id,
            design_run_id = design$run_id,
            import_run_id = import$run_id
        )
    }
    invisible(TRUE)
}

protDiaArtifactNormalizedRefs <- function(refs) {
    normalized <- lapply(refs, artifactStoreNormalizeRef)
    names(normalized) <- names(refs)
    normalized
}

protDiaArtifactValidateStateEvidence <- function(manager, evidence, state_object) {
    metadata <- manager$getStateMetadata()
    audit <- metadata$audit_metadata
    expected_refs <- protDiaArtifactNormalizedRefs(evidence$design$refs)
    actual_refs <- protDiaArtifactNormalizedRefs(audit$stage_artifact_refs)
    valid <- identical(
        manager$getCurrentGenerationId(),
        evidence$current_state$generation_id
    ) && identical(
        manager$getCurrentStateName(),
        evidence$design$parameters$state_name
    ) && identical(audit$stage_id, "design") &&
        identical(audit$run_id, evidence$design$run_id) &&
        protDiaArtifactContractVersion(audit$readthrough_contract_version) &&
        identical(audit$parent_import_run_id, evidence$import$run_id) &&
        identical(
            audit$parent_import_generation_id,
            evidence$import$generation_id
        ) && identical(
            audit$parent_import_artifact_id,
            evidence$import$refs$canonical_data$artifact_id
        ) && identical(
            audit$parent_import_semantic_digest,
            evidence$import$refs$canonical_data$hash_policy$semantic$digest
        ) && identical(actual_refs, expected_refs) && isS4(state_object) &&
        identical(class(state_object)[[1L]], "PeptideQuantitativeData") &&
        identical(methods::validObject(state_object, test = TRUE), TRUE)
    if (!isTRUE(valid)) {
        protDiaArtifactAbort(
            "DIA-NN current S4 state does not match its committed design run",
            "multischolar_prot_dia_state_lineage_mismatch",
            design_run_id = evidence$design$run_id,
            state_generation_id = evidence$current_state$generation_id
        )
    }
    invisible(TRUE)
}

protDiaArtifactRestoreOptional <- function(value, available, role) {
    if (!is.logical(available) || length(available) != 1L || is.na(available)) {
        protDiaArtifactAbort(
            sprintf("DIA-NN '%s' availability flag is invalid", role),
            "multischolar_invalid_prot_dia_run_parameters"
        )
    }
    if (isTRUE(available)) return(value)
    if (!identical(value, protDiaArtifactOptionalTable(NULL, role))) {
        protDiaArtifactAbort(
            sprintf("DIA-NN '%s' absence marker is invalid", role),
            "multischolar_inexact_prot_dia_artifact_hydration"
        )
    }
    NULL
}

protDiaArtifactRestoreContrasts <- function(value, kind) {
    if (identical(kind, "data.frame")) return(value)
    if (identical(kind, "character") && identical(names(value), "contrasts")) {
        return(as.character(value$contrasts))
    }
    if (identical(kind, "null") &&
        identical(value, protDiaArtifactOptionalTable(NULL, "contrasts"))) {
        return(NULL)
    }
    protDiaArtifactAbort(
        "DIA-NN contrast representation is incompatible with its provenance",
        "multischolar_inexact_prot_dia_artifact_hydration"
    )
}

protDiaArtifactValidateScientificTables <- function(
    state_object,
    config,
    design_tables
) {
    expected_args <- protDiaArtifactMetadataTable(state_object@args, "args")
    valid <- identical(design_tables$cleaned_data, state_object@peptide_data) &&
        identical(design_tables$design_matrix, state_object@design_matrix) &&
        identical(design_tables$args, expected_args) &&
        identical(config, state_object@args)
    if (!isTRUE(valid)) {
        protDiaArtifactAbort(
            "DIA-NN design tables differ from the current scientific S4 state",
            "multischolar_inexact_prot_dia_artifact_hydration"
        )
    }
    invisible(TRUE)
}

protDiaArtifactValidateSettledScientificTables <- function(
    state_object,
    config,
    design_tables
) {
    expected_args <- protDiaArtifactMetadataTable(state_object@args, "args")
    valid <- identical(design_tables$design_matrix, state_object@design_matrix) &&
        identical(design_tables$args, expected_args) &&
        identical(config, state_object@args)
    if (!isTRUE(valid)) {
        protDiaArtifactAbort(
            "DIA-NN settled metadata differs from the current scientific S4 state",
            "multischolar_inexact_prot_dia_artifact_hydration"
        )
    }
    invisible(TRUE)
}

protDiaArtifactColumnMapping <- function(parameters) {
    encoded <- parameters$column_mapping_serialized
    mapping <- artifactWorkflowStateUnserializeMetadata(
        encoded,
        "DIA-NN import column mapping"
    )
    if (!is.list(mapping)) {
        protDiaArtifactAbort(
            "DIA-NN import column mapping is not an R list",
            "multischolar_invalid_prot_dia_run_parameters"
        )
    }
    mapping
}

newProtDiaArtifactStateManager <- function(context, resource_policy = NULL) {
    newWorkflowState(
        workflow_context = context,
        resource_policy = resource_policy,
        workflow_descriptor = artifactDiaWorkflowDescriptor(),
        descriptor_catalogue = artifactWorkflowDescriptorCatalogue(),
        codec_catalogue = artifactS4CodecCatalogue()
    )
}

hydrateProtDiaResumeBundle <- function(
    context,
    resource_policy = NULL,
    retain_source_payloads = TRUE
) {
    if (!is.logical(retain_source_payloads) ||
        length(retain_source_payloads) != 1L || is.na(retain_source_payloads)) {
        protDiaArtifactAbort(
            "DIA-NN source retention mode must be true or false",
            "multischolar_invalid_prot_dia_readthrough_mode"
        )
    }
    evidence <- collectProtDiaResumeEvidence(
        context,
        resource_policy,
        payload_validation = if (isTRUE(retain_source_payloads)) {
            "sidecar"
        } else {
            "digest"
        }
    )
    protDiaArtifactValidateResumeContract(evidence)
    if (isTRUE(retain_source_payloads)) {
        import_tables <- readProtDiaArtifactStage(evidence$store, evidence$import)
        design_tables <- readProtDiaArtifactStage(evidence$store, evidence$design)
    } else {
        import_tables <- list(canonical_data = NULL)
        design_tables <- readProtDiaArtifactStageRoles(
            evidence$store,
            evidence$design,
            c("design_matrix", "contrasts", "args", "annotations", "sequences")
        )
        design_tables$cleaned_data <- NULL
    }
    manager <- newProtDiaArtifactStateManager(context, resource_policy)
    manager_owned <- TRUE
    on.exit({
        if (manager_owned) try(manager$close(), silent = TRUE)
    }, add = TRUE)
    state_object <- manager$getState()
    protDiaArtifactValidateStateEvidence(manager, evidence, state_object)
    config <- manager$getStateConfig()
    if (isTRUE(retain_source_payloads)) {
        protDiaArtifactValidateScientificTables(state_object, config, design_tables)
    } else {
        protDiaArtifactValidateSettledScientificTables(
            state_object,
            config,
            design_tables
        )
    }
    parameters <- evidence$design$parameters
    bundle <- list(
        context = context,
        state_manager = manager,
        state_object = state_object,
        data_tbl = import_tables$canonical_data,
        data_cln = design_tables$cleaned_data,
        design_matrix = design_tables$design_matrix,
        contrasts_tbl = protDiaArtifactRestoreContrasts(
            design_tables$contrasts,
            parameters$contrasts_kind
        ),
        config_list = config,
        column_mapping = protDiaArtifactColumnMapping(evidence$import$parameters),
        uniprot_dat_cln = protDiaArtifactRestoreOptional(
            design_tables$annotations,
            parameters$annotation_available,
            "annotations"
        ),
        aa_seq_tbl_final = protDiaArtifactRestoreOptional(
            design_tables$sequences,
            parameters$sequence_data_available,
            "sequences"
        ),
        source_payloads_retained = retain_source_payloads,
        annotation_completed = TRUE,
        readthrough_mode = if (isTRUE(retain_source_payloads)) {
            "full"
        } else {
            "settled"
        },
        evidence = evidence
    )
    manager_owned <- FALSE
    bundle
}

applyProtDiaResumeBundle <- function(workflow_data, bundle) {
    previous_manager <- workflow_data$state_manager
    workflow_data$workflow_context <- bundle$context
    workflow_data$state_manager <- bundle$state_manager
    workflow_data$data_tbl <- bundle$data_tbl
    workflow_data$data_cln <- bundle$data_cln
    workflow_data$design_matrix <- bundle$design_matrix
    workflow_data$contrasts_tbl <- bundle$contrasts_tbl
    workflow_data$config_list <- bundle$config_list
    workflow_data$column_mapping <- bundle$column_mapping
    workflow_data$uniprot_dat_cln <- bundle$uniprot_dat_cln
    workflow_data$aa_seq_tbl_final <- bundle$aa_seq_tbl_final
    workflow_data$data_format <- "diann"
    workflow_data$data_type <- "peptide"
    workflow_data$artifact_readthrough_refs <- list(
        import = bundle$evidence$import$refs,
        design = bundle$evidence$design$refs
    )
    readthrough <- recordProtDiaReadthroughProof(workflow_data, bundle)
    status <- workflow_data$tab_status
    status$setup_import <- "complete"
    status$design_matrix <- "complete"
    if (identical(status$quality_control, "disabled")) {
        status$quality_control <- "pending"
    }
    workflow_data$tab_status <- status
    if (inherits(previous_manager, "ArtifactWorkflowState") &&
        !identical(previous_manager, bundle$state_manager)) {
        try(previous_manager$close(), silent = TRUE)
    }
    invisible(readthrough)
}

protDiaArtifactLegacyEligibleError <- function(error) {
    inherits(error, c(
        "multischolar_missing_prot_dia_committed_stage",
        "multischolar_incomplete_prot_dia_readthrough_contract",
        "multischolar_missing_prot_dia_current_state"
    ))
}

resumeProtDiaArtifactWorkflow <- function(
    workflow_data,
    experiment_paths,
    experiment_label,
    storage_policy = NULL,
    resource_policy = NULL,
    legacy_adapter = NULL
) {
    prepared <- createProtDiaResumeContext(
        experiment_paths,
        experiment_label,
        storage_policy
    )
    if (!isTRUE(prepared$enabled)) {
        return(c(prepared, list(ok = TRUE, resumed = FALSE)))
    }
    retain_source_payloads <- !identical(
        prepared$context$getStorageDecision()$effective_rollout,
        "evict"
    )
    bundle <- tryCatch(
        hydrateProtDiaResumeBundle(
            prepared$context,
            resource_policy,
            retain_source_payloads = retain_source_payloads
        ),
        error = \(error) error
    )
    if (inherits(bundle, "error")) {
        if (protDiaArtifactLegacyEligibleError(bundle) &&
            is.function(legacy_adapter)) {
            legacy <- legacy_adapter(
                workflow_data = workflow_data,
                experiment_paths = experiment_paths,
                error = bundle
            )
            return(list(
                enabled = FALSE,
                ok = TRUE,
                resumed = FALSE,
                reason = "explicit_legacy_adapter",
                legacy_result = legacy
            ))
        }
        stop(bundle)
    }
    applied <- FALSE
    on.exit({
        if (!applied) try(bundle$state_manager$close(), silent = TRUE)
    }, add = TRUE)
    readthrough <- applyProtDiaResumeBundle(workflow_data, bundle)
    applied <- TRUE
    list(
        enabled = TRUE,
        ok = TRUE,
        resumed = TRUE,
        reason = "artifact_readthrough",
        project_id = bundle$evidence$identity$project_id,
        import_run_id = bundle$evidence$import$run_id,
        design_run_id = bundle$evidence$design$run_id,
        state_generation_id = bundle$evidence$current_state$generation_id,
        import_ref = unclass(bundle$evidence$import$refs$canonical_data),
        source_payloads_retained = bundle$source_payloads_retained,
        readthrough_mode = bundle$readthrough_mode,
        compatibility_checkpoint = readthrough$compatibility_checkpoint
    )
}

protDiaArtifactProjectPresent <- function(
    experiment_paths,
    experiment_label,
    storage_policy
) {
    policy <- normalizeWorkflowStoragePolicy(storage_policy)
    if (identical(policy$requested_backend, "memory")) return(FALSE)
    tryCatch(
        isTRUE(protDiaArtifactProbeManifest(
            experiment_paths,
            experiment_label
        )$found),
        error = \(...) TRUE
    )
}

resumeProtDiaArtifactWorkflowSafely <- function(
    workflow_data,
    experiment_paths,
    experiment_label,
    storage_policy = NULL,
    resource_policy = NULL,
    log_warn = logger::log_warn
) {
    artifact_project <- protDiaArtifactProjectPresent(
        experiment_paths,
        experiment_label,
        storage_policy
    )
    result <- tryCatch(
        resumeProtDiaArtifactWorkflow(
            workflow_data,
            experiment_paths,
            experiment_label,
            storage_policy,
            resource_policy
        ),
        error = \(error) {
            log_warn(paste(
                "DIA-NN artifact read-through failed without changing workflow state:",
                conditionMessage(error)
            ))
            list(
                enabled = artifact_project,
                ok = FALSE,
                resumed = FALSE,
                reason = "artifact_readthrough_failed",
                error_class = class(error),
                error_message = conditionMessage(error),
                remediation = error$remediation
            )
        }
    )
    result$artifact_project <- artifact_project
    recordProtDiaArtifactResult(workflow_data, "resume", result)
    result
}

previewProtDiaImportArtifact <- function(
    context,
    import_ref = NULL,
    projections = NULL,
    filters = list(),
    limit = NULL,
    resource_policy = NULL,
    query_session = NULL
) {
    if (is.null(import_ref)) {
        evidence <- collectProtDiaResumeEvidence(context, resource_policy)
        protDiaArtifactValidateResumeContract(evidence)
        store <- evidence$store
        import_ref <- evidence$import$refs$canonical_data
    } else {
        if (!inherits(context, "WorkflowContext") || !context$isBound() ||
            !protDiaArtifactIdentityMatches(context$getIdentity())) {
            protDiaArtifactAbort(
                "DIA-NN preview requires its exact bound artifact context",
                "multischolar_invalid_prot_dia_artifact_context"
            )
        }
        identity <- context$getIdentity()
        store <- newArtifactStore(context$getPaths(), identity$project_id)
    }
    descriptor <- artifactDiaWorkflowDescriptor()
    queryArtifactRef(
        store,
        import_ref,
        descriptor,
        names(descriptor$queries)[[1L]],
        projections = projections,
        filters = filters,
        limit = limit,
        resource_policy = resource_policy,
        query_session = query_session
    )
}
