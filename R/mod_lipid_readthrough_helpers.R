# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

.LIPID_READTHROUGH_PROOF_SCHEMA <- "multischolar.lipidomics.readthrough_proof"
.LIPID_READTHROUGH_PROOF_VERSION <- 1L

lipidArtifactContractVersion <- function(value) {
    artifactStageReadthroughContractVersion(
        value,
        .LIPID_READTHROUGH_VERSION
    )
}

validateLipidResumeContract <- function(evidence) {
    adapter <- lipidArtifactReadthroughAdapter()
    import <- evidence$import
    design <- evidence$design
    import_parameters <- import$parameters
    design_parameters <- design$parameters
    import_order <- unlist(import_parameters$assay_order, use.names = FALSE)
    design_order <- unlist(design_parameters$assay_order, use.names = FALSE)
    valid <- lipidArtifactContractVersion(
        import_parameters$readthrough_contract_version
    ) && lipidArtifactContractVersion(
        design_parameters$readthrough_contract_version
    ) && identical(
        import_parameters$capability_id,
        evidence$descriptor$descriptor_id
    ) && identical(
        design_parameters$capability_id,
        evidence$descriptor$descriptor_id
    ) && identical(import_parameters$input_format, "lipidsearch") &&
        identical(import_parameters$data_level, "lipid") &&
        identical(design_parameters$state_name, "lipid_raw_data_s4") &&
        identical(design_parameters$workflow_type, "lipidomics_standard") &&
        identical(import_order, design_order) &&
        identical(design_parameters$parent_import_run_id, import$run_id) &&
        identical(
            design_parameters$parent_import_generation_id,
            import$generation_id
        ) && identical(
            evidence$current_state$logical_name,
            design_parameters$state_name
        ) && artifactStageIdentityMatches(
            evidence$identity,
            evidence$descriptor
        )
    if (!isTRUE(valid)) {
        artifactStageReadthroughAbort(
            adapter,
            "incomplete_contract",
            "lipidomics import, design, and state are not one exact lineage"
        )
    }
    invisible(TRUE)
}

newLipidResumeStateManager <- function(
    context,
    resource_policy = NULL,
    settled_bootstrap = NULL
) {
    ArtifactWorkflowState$new(
        workflow_context = context,
        resource_policy = resource_policy,
        dehydrate_fn = dehydrateLipidomicsS4Artifact,
        validate_bundle_fn = validateLipidomicsS4Bundle,
        hydrate_fn = hydrateLipidomicsS4Artifact,
        descriptor_contract = artifactStageDescriptorContract(
            artifactLipidomicsWorkflowDescriptor()
        ),
        settled_bootstrap = settled_bootstrap
    )
}

validateLipidStateEvidence <- function(manager, evidence, state_object) {
    adapter <- lipidArtifactReadthroughAdapter()
    metadata <- manager$getStateMetadata()
    audit <- metadata$audit_metadata
    expected_refs <- artifactPayloadReadthroughRefProof(evidence$design$refs)
    actual_refs <- artifactPayloadReadthroughRefProof(audit$stage_artifact_refs)
    valid <- identical(
        manager$getCurrentGenerationId(),
        evidence$current_state$generation_id
    ) && identical(
        manager$getCurrentStateName(),
        evidence$design$parameters$state_name
    ) && identical(audit$stage_id, "design") &&
        identical(audit$run_id, evidence$design$run_id) &&
        lipidArtifactContractVersion(audit$readthrough_contract_version) &&
        identical(audit$parent_import_run_id, evidence$import$run_id) &&
        identical(
            audit$parent_import_generation_id,
            evidence$import$generation_id
        ) && identical(actual_refs, expected_refs) &&
        methods::is(state_object, "LipidomicsAssayData") &&
        identical(methods::validObject(state_object, test = TRUE), TRUE)
    if (!isTRUE(valid)) {
        artifactStageReadthroughAbort(
            adapter,
            "incomplete_contract",
            "lipidomics current S4 state differs from its design generation"
        )
    }
    invisible(TRUE)
}

validateLipidSettledStateEvidence <- function(manager, evidence) {
    adapter <- lipidArtifactReadthroughAdapter()
    metadata <- manager$getStateMetadata()
    audit <- metadata$audit_metadata
    expected_refs <- artifactPayloadReadthroughRefProof(evidence$design$refs)
    actual_refs <- artifactPayloadReadthroughRefProof(audit$stage_artifact_refs)
    valid <- identical(
        manager$getCurrentGenerationId(),
        evidence$current_state$generation_id
    ) && identical(
        manager$getCurrentStateName(),
        evidence$design$parameters$state_name
    ) && identical(metadata$s4_class, "LipidomicsAssayData") &&
        identical(audit$stage_id, "design") &&
        identical(audit$run_id, evidence$design$run_id) &&
        lipidArtifactContractVersion(audit$readthrough_contract_version) &&
        identical(audit$parent_import_run_id, evidence$import$run_id) &&
        identical(
            audit$parent_import_generation_id,
            evidence$import$generation_id
        ) && identical(actual_refs, expected_refs)
    if (!isTRUE(valid)) {
        artifactStageReadthroughAbort(
            adapter,
            "incomplete_contract",
            "lipidomics settled state metadata differs from its design generation"
        )
    }
    invisible(TRUE)
}

lipidResumeAssayTables <- function(tables, parameters, prefix = "assay") {
    assay_names <- unlist(parameters$assay_order, use.names = FALSE)
    roles <- if (prefix == "assay") {
        unlist(parameters$assay_roles, use.names = FALSE)
    } else {
        sprintf("cleaned_assay_%04d", seq_along(assay_names))
    }
    values <- tables[roles]
    names(values) <- assay_names
    values
}

validateLipidImportTables <- function(evidence, tables, mapping) {
    adapter <- lipidArtifactReadthroughAdapter()
    parameters <- evidence$import$parameters
    assays <- lipidResumeAssayTables(tables, parameters)
    manifest <- tables$assay_manifest
    assay_names <- names(assays)
    roles <- unlist(parameters$assay_roles, use.names = FALSE)
    valid <- is.data.frame(manifest) &&
        identical(manifest$assay_name, assay_names) &&
        identical(manifest$assay_order, seq_along(assay_names)) &&
        identical(manifest$artifact_role, roles) &&
        identical(
            tables$column_mapping,
            artifactStageMetadataTable(mapping)
        ) && identical(
            tables$source_manifest$assay_name,
            assay_names
        )
    if (isTRUE(valid)) {
        valid <- all(vapply(seq_along(assays), \(index) {
            assay <- assays[[index]]
            identical(
                artifactSemanticDigest(assay),
                manifest$table_digest[[index]]
            ) && identical(nrow(assay), manifest$feature_count[[index]]) &&
                identical(
                    length(mapping$sample_columns),
                    manifest$sample_count[[index]]
                ) && identical(
                    names(assay),
                    jsonlite::fromJSON(manifest$column_schema_json[[index]])
                )
        }, logical(1)))
    }
    if (!isTRUE(valid)) {
        artifactStageReadthroughAbort(
            adapter,
            "incomplete_contract",
            "lipidomics import tables differ from their assay manifest"
        )
    }
    assays
}

validateLipidDesignTables <- function(
    evidence,
    tables,
    mapping,
    state_object,
    config
) {
    adapter <- lipidArtifactReadthroughAdapter()
    parameters <- evidence$design$parameters
    assays <- lipidResumeAssayTables(tables, parameters, "cleaned")
    alignment <- tables$assay_alignment
    dependencies <- tables$raw_s4_dependencies
    state_assays <- methods::slot(state_object, "lipid_data")
    valid <- identical(assays, state_assays) && identical(
        tables$design_matrix,
        methods::slot(state_object, "design_matrix")
    ) && identical(
        tables$args,
        artifactStageMetadataTable(methods::slot(state_object, "args"))
    ) && identical(config, methods::slot(state_object, "args")) &&
        identical(tables$column_mapping, artifactStageMetadataTable(mapping)) &&
        identical(alignment$assay_name, names(state_assays)) &&
        identical(alignment$assay_order, seq_along(state_assays)) &&
        identical(dependencies$assay_name, names(state_assays)) &&
        all(dependencies$parent_import_run_id == evidence$import$run_id) &&
        all(dependencies$parent_import_generation_id == evidence$import$generation_id)
    if (isTRUE(valid)) {
        valid <- all(vapply(seq_along(state_assays), \(index) {
            assay <- state_assays[[index]]
            feature_column <- methods::slot(state_object, "lipid_id_column")
            identical(
                alignment$feature_order_digest[[index]],
                artifactSemanticDigest(assay[[feature_column]])
            )
        }, logical(1)))
    }
    if (!isTRUE(valid)) {
        artifactStageReadthroughAbort(
            adapter,
            "incomplete_contract",
            "lipidomics design tables differ from the current S4 state"
        )
    }
    assays
}

validateLipidSettledTables <- function(
    evidence,
    import_tables,
    design_tables,
    mapping,
    config
) {
    adapter <- lipidArtifactReadthroughAdapter()
    assay_order <- unlist(
        evidence$design$parameters$assay_order,
        use.names = FALSE
    )
    assay_roles <- unlist(
        evidence$import$parameters$assay_roles,
        use.names = FALSE
    )
    manifest <- import_tables$assay_manifest
    alignment <- design_tables$assay_alignment
    dependencies <- design_tables$raw_s4_dependencies
    sample_columns <- mapping$sample_columns
    valid <- is.data.frame(manifest) && is.data.frame(alignment) &&
        is.data.frame(dependencies) && is.data.frame(design_tables$design_matrix) &&
        identical(manifest$assay_name, assay_order) &&
        identical(manifest$assay_order, seq_along(assay_order)) &&
        identical(manifest$artifact_role, assay_roles) &&
        all(manifest$sample_count == length(sample_columns)) &&
        identical(import_tables$source_manifest$assay_name, assay_order) &&
        identical(
            import_tables$column_mapping,
            artifactStageMetadataTable(mapping)
        ) && identical(
            design_tables$column_mapping,
            artifactStageMetadataTable(mapping)
        ) && identical(design_tables$args, artifactStageMetadataTable(config)) &&
        identical(alignment$assay_name, assay_order) &&
        identical(alignment$assay_order, seq_along(assay_order)) &&
        identical(dependencies$assay_name, assay_order) &&
        all(dependencies$parent_import_run_id == evidence$import$run_id) &&
        all(
            dependencies$parent_import_generation_id ==
                evidence$import$generation_id
        ) && identical(
            as.character(design_tables$design_matrix$Run),
            as.character(sample_columns)
        )
    if (!isTRUE(valid)) {
        artifactStageReadthroughAbort(
            adapter,
            "incomplete_contract",
            "lipidomics settled metadata differs from committed assay lineage"
        )
    }
    invisible(TRUE)
}

hydrateLipidResumeBundle <- function(
    context,
    resource_policy = NULL,
    retain_source_payloads = TRUE
) {
    adapter <- lipidArtifactReadthroughAdapter()
    if (!is.logical(retain_source_payloads) ||
        length(retain_source_payloads) != 1L || is.na(retain_source_payloads)) {
        artifactStageReadthroughAbort(
            adapter,
            "incomplete_contract",
            "lipidomics source retention mode must be true or false"
        )
    }
    if (!isTRUE(retain_source_payloads)) {
        resource_policy <- artifactPayloadSettledRegistryPolicy(resource_policy)
    }
    snapshot <- if (isTRUE(retain_source_payloads)) {
        NULL
    } else {
        artifactSettledResumeRead(
            context,
            artifactLipidomicsWorkflowDescriptor(),
            adapter
        )
    }
    if (!is.null(snapshot)) artifactReleaseTransientMemory()
    evidence <- if (is.null(snapshot)) {
        collectLipidResumeEvidence(
            context,
            resource_policy,
            payload_validation = if (isTRUE(retain_source_payloads)) {
                "sidecar"
            } else {
                "digest"
            }
        )
    } else {
        list(
            identity = snapshot$identity,
            descriptor = artifactLipidomicsWorkflowDescriptor(),
            store = snapshot$store,
            import = snapshot$import,
            design = snapshot$design,
            current_state = snapshot$current_state
        )
    }
    validateLipidResumeContract(evidence)
    if (isTRUE(retain_source_payloads)) {
        import_tables <- readArtifactStage(
            adapter,
            evidence$store,
            evidence$import
        )
        design_tables <- readArtifactStage(
            adapter,
            evidence$store,
            evidence$design
        )
    } else if (!is.null(snapshot)) {
        import_tables <- snapshot$settled_tables$import
        design_tables <- snapshot$settled_tables$design
    } else {
        import_tables <- readArtifactStageRoles(
            adapter,
            evidence$store,
            evidence$import,
            .LIPID_IMPORT_FIXED_ROLES
        )
        design_tables <- readArtifactStageRoles(
            adapter,
            evidence$store,
            evidence$design,
            .LIPID_DESIGN_FIXED_ROLES
        )
    }
    manager <- newLipidResumeStateManager(
        context,
        resource_policy,
        settled_bootstrap = if (is.null(snapshot)) {
            NULL
        } else {
            snapshot$state_bootstrap
        }
    )
    manager_owned <- TRUE
    on.exit({
        if (manager_owned) try(manager$close(), silent = TRUE)
    }, add = TRUE)
    config <- manager$getStateConfig()
    mapping <- artifactStageColumnMapping(
        adapter,
        evidence$import$parameters
    )
    if (isTRUE(retain_source_payloads)) {
        state_object <- manager$getState()
        validateLipidStateEvidence(manager, evidence, state_object)
        data_tbl <- validateLipidImportTables(evidence, import_tables, mapping)
        data_cln <- validateLipidDesignTables(
            evidence,
            design_tables,
            mapping,
            state_object,
            config
        )
    } else {
        state_object <- NULL
        validateLipidSettledStateEvidence(manager, evidence)
        validateLipidSettledTables(
            evidence,
            import_tables,
            design_tables,
            mapping,
            config
        )
        data_tbl <- NULL
        data_cln <- NULL
    }
    contrasts <- artifactStageRestoreContrasts(
        adapter,
        design_tables$contrasts,
        evidence$design$parameters$contrasts_kind
    )
    if (isTRUE(retain_source_payloads)) {
        preflight <- validateLipidDesignDaPreflight(
            designMatrix = design_tables$design_matrix,
            assayList = data_cln,
            contrastsTbl = contrasts,
            formulaString = evidence$design$parameters$formula_string,
            columnMapping = mapping,
            requireContrasts = FALSE
        )
        if (!isTRUE(preflight$valid)) {
            artifactStageReadthroughAbort(
                adapter,
                "incomplete_contract",
                "lipidomics resumed design fails scientific preflight"
            )
        }
    }
    if (!isTRUE(retain_source_payloads)) artifactReleaseTransientMemory()
    bundle <- list(
        context = context,
        descriptor = evidence$descriptor,
        state_manager = manager,
        state_object = state_object,
        data_tbl = data_tbl,
        data_cln = data_cln,
        design_matrix = design_tables$design_matrix,
        contrasts_tbl = contrasts,
        config_list = config,
        column_mapping = mapping,
        source_payloads_retained = retain_source_payloads,
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

lipidReadthroughConsumerInventory <- function() {
    readers <- c(
        "previewLipidImportArtifact", "validateLipidDesignDaPreflight",
        "runLipidQcS4FinalizeWorkflow", "normaliseBetweenSamples",
        "handleLipidNormExportSession", "differentialAbundanceAnalysisHelper",
        "getLipidSummaryStateObject"
    )
    data.frame(
        category = c(
            "preview", "compatibility", "qc", "normalization", "session",
            "da", "report"
        ),
        reader = readers,
        restored_source = c(
            "ordered immutable import refs", "design and contrast tables",
            "current exact S4 state", "current exact S4 state",
            "state plus coordinator metadata", "assay/design/contrast state",
            "current state and workflow metadata"
        ),
        verified = vapply(readers, \(reader) exists(
            reader,
            envir = environment(lipidReadthroughConsumerInventory),
            mode = "function",
            inherits = TRUE
        ), logical(1)),
        stringsAsFactors = FALSE
    )
}

validateLipidReadthroughConsumerInventory <- function(inventory) {
    categories <- c(
        "preview", "compatibility", "qc", "normalization", "session", "da",
        "report"
    )
    valid <- is.data.frame(inventory) && identical(
        names(inventory),
        c("category", "reader", "restored_source", "verified")
    ) && setequal(unique(inventory$category), categories) &&
        !anyDuplicated(inventory$reader) && all(inventory$verified)
    if (!isTRUE(valid)) {
        lipidArtifactAbort(
            "lipidomics read-through consumer inventory is incomplete",
            "multischolar_incomplete_lipidomics_consumer_inventory"
        )
    }
    inventory
}

lipidReadthroughResourceSnapshot <- function(workflow_data, phase) {
    fields <- c("data_tbl", "data_cln", "design_matrix")
    object_bytes <- vapply(fields, \(field) {
        value <- tryCatch(workflow_data[[field]], error = \(...) NULL)
        unname(as.numeric(utils::object.size(value)))
    }, numeric(1))
    manager <- tryCatch(workflow_data$state_manager, error = \(...) NULL)
    cache <- tryCatch(manager$getCacheInfo(), error = \(...) NULL)
    state <- if (is.list(cache) && cache$entries > 0L) {
        tryCatch(manager$getState(), error = \(...) NULL)
    } else if (inherits(manager, "WorkflowState")) {
        tryCatch(manager$getState(), error = \(...) NULL)
    } else {
        NULL
    }
    list(
        phase = phase,
        rss_bytes = unname(as.numeric(
            ps::ps_memory_info(ps::ps_handle())[["rss"]]
        )),
        object_bytes = object_bytes,
        state_object_bytes = unname(as.numeric(utils::object.size(state))),
        captured_at = artifactRefUtcNow()
    )
}

newLipidReadthroughProof <- function(bundle, inventory) {
    evidence <- bundle$evidence
    list(
        schema = .LIPID_READTHROUGH_PROOF_SCHEMA,
        schema_version = .LIPID_READTHROUGH_PROOF_VERSION,
        descriptor_id = evidence$descriptor$descriptor_id,
        descriptor_version = evidence$descriptor$descriptor_version,
        project_id = evidence$identity$project_id,
        workflow_id = evidence$identity$workflow_id,
        import_run_id = evidence$import$run_id,
        design_run_id = evidence$design$run_id,
        state_generation_id = evidence$current_state$generation_id,
        assay_order = unlist(
            evidence$design$parameters$assay_order,
            use.names = FALSE
        ),
        payload_validated = TRUE,
        registry_ready = TRUE,
        current_pointer_verified = TRUE,
        readthrough_completed = TRUE,
        annotation_completed = TRUE,
        readthrough_mode = bundle$readthrough_mode,
        source_payloads_retained = bundle$source_payloads_retained,
        eviction_performed = !bundle$source_payloads_retained,
        import_refs = artifactPayloadReadthroughRefProof(evidence$import$refs),
        design_refs = artifactPayloadReadthroughRefProof(evidence$design$refs),
        consumer_inventory_digest = artifactSemanticDigest(inventory),
        verified_at = artifactRefUtcNow()
    )
}

captureLipidResumeFields <- function(workflow_data, fields) {
    values <- lapply(fields, \(field) {
        tryCatch(workflow_data[[field]], error = \(...) NULL)
    })
    stats::setNames(values, fields)
}

restoreLipidResumeFields <- function(workflow_data, values) {
    for (field in names(values)) workflow_data[[field]] <- values[[field]]
    invisible(TRUE)
}

applyLipidResumeBundle <- function(
    workflow_data,
    bundle,
    failure_injector = NULL,
    inventory_fn = lipidReadthroughConsumerInventory
) {
    inventory <- validateLipidReadthroughConsumerInventory(inventory_fn())
    fields <- c(
        "workflow_context", "state_manager", "data_tbl", "data_cln",
        "design_matrix", "contrasts_tbl", "config_list", "column_mapping",
        "data_format", "data_type", "artifact_readthrough_refs",
        "artifact_readthrough_proof", "artifact_compatibility_checkpoint",
        "tab_status"
    )
    previous <- captureLipidResumeFields(workflow_data, fields)
    committed <- FALSE
    on.exit({
        if (!committed) restoreLipidResumeFields(workflow_data, previous)
    }, add = TRUE)
    artifactStoreInvokeFailure(
        failure_injector,
        "before_resume_apply",
        list(descriptor_id = bundle$descriptor$descriptor_id)
    )
    workflow_data$workflow_context <- bundle$context
    workflow_data$state_manager <- bundle$state_manager
    workflow_data$data_tbl <- bundle$data_tbl
    workflow_data$data_cln <- bundle$data_cln
    workflow_data$design_matrix <- bundle$design_matrix
    workflow_data$contrasts_tbl <- bundle$contrasts_tbl
    workflow_data$config_list <- bundle$config_list
    workflow_data$column_mapping <- bundle$column_mapping
    workflow_data$data_format <- "lipidsearch"
    workflow_data$data_type <- "lipid"
    workflow_data$artifact_readthrough_refs <- list(
        import = bundle$evidence$import$refs,
        design = bundle$evidence$design$refs
    )
    proof <- newLipidReadthroughProof(bundle, inventory)
    checkpoint <- list(
        strategy = "restore_memory_reads_keep_lipidomics_generations",
        verified = TRUE,
        descriptor_id = proof$descriptor_id,
        project_id = proof$project_id,
        workflow_id = proof$workflow_id,
        design_run_id = proof$design_run_id,
        state_generation_id = proof$state_generation_id,
        consumer_inventory_digest = proof$consumer_inventory_digest,
        verified_at = proof$verified_at
    )
    workflow_data$artifact_readthrough_proof <- proof
    workflow_data$artifact_compatibility_checkpoint <- checkpoint
    status <- workflow_data$tab_status
    status$setup_import <- "complete"
    status$design_matrix <- "complete"
    if (identical(status$quality_control, "disabled")) {
        status$quality_control <- "pending"
    }
    workflow_data$tab_status <- status
    artifactStoreInvokeFailure(
        failure_injector,
        "after_resume_apply",
        list(descriptor_id = bundle$descriptor$descriptor_id)
    )
    committed <- TRUE
    previous_manager <- previous$state_manager
    if (inherits(previous_manager, "ArtifactWorkflowState") &&
        !identical(previous_manager, bundle$state_manager)) {
        try(previous_manager$close(), silent = TRUE)
    }
    invisible(list(proof = proof, compatibility_checkpoint = checkpoint))
}

resumeLipidArtifactWorkflow <- function(
    workflow_data,
    experiment_paths,
    experiment_label,
    storage_policy = NULL,
    resource_policy = NULL,
    failure_injector = NULL,
    inventory_fn = lipidReadthroughConsumerInventory
) {
    if (is.null(storage_policy) &&
        inherits(workflow_data$workflow_context, "WorkflowContext")) {
        storage_policy <- workflow_data$workflow_context$getStoragePolicy()
    }
    prepared <- createLipidResumeContext(
        experiment_paths,
        experiment_label,
        storage_policy
    )
    if (!isTRUE(prepared$enabled)) {
        return(c(prepared, list(ok = TRUE, resumed = FALSE)))
    }
    before <- lipidReadthroughResourceSnapshot(workflow_data, "before")
    retain_source_payloads <- !identical(
        prepared$context$getStorageDecision()$effective_rollout,
        "evict"
    )
    bundle <- hydrateLipidResumeBundle(
        prepared$context,
        resource_policy,
        retain_source_payloads = retain_source_payloads
    )
    applied <- FALSE
    on.exit({
        if (!applied) try(bundle$state_manager$close(), silent = TRUE)
    }, add = TRUE)
    readthrough <- applyLipidResumeBundle(
        workflow_data,
        bundle,
        failure_injector,
        inventory_fn
    )
    applied <- TRUE
    after <- lipidReadthroughResourceSnapshot(workflow_data, "after")
    list(
        enabled = TRUE,
        ok = TRUE,
        resumed = TRUE,
        reason = "artifact_readthrough",
        capability_id = bundle$descriptor$descriptor_id,
        project_id = bundle$evidence$identity$project_id,
        import_run_id = bundle$evidence$import$run_id,
        design_run_id = bundle$evidence$design$run_id,
        state_generation_id = bundle$evidence$current_state$generation_id,
        assay_order = unlist(
            bundle$evidence$design$parameters$assay_order,
            use.names = FALSE
        ),
        source_payloads_retained = bundle$source_payloads_retained,
        readthrough_mode = bundle$readthrough_mode,
        eviction_performed = !bundle$source_payloads_retained,
        compatibility_checkpoint = readthrough$compatibility_checkpoint,
        resource_evidence = list(before = before, after = after)
    )
}

resumeLipidArtifactWorkflowSafely <- function(
    workflow_data,
    experiment_paths,
    experiment_label,
    storage_policy = NULL,
    resource_policy = NULL,
    failure_injector = NULL,
    inventory_fn = lipidReadthroughConsumerInventory,
    log_warn = logger::log_warn
) {
    artifact_project <- lipidArtifactProjectPresent(
        experiment_paths,
        experiment_label,
        storage_policy
    )
    result <- tryCatch(
        resumeLipidArtifactWorkflow(
            workflow_data,
            experiment_paths,
            experiment_label,
            storage_policy,
            resource_policy,
            failure_injector,
            inventory_fn
        ),
        error = \(error) {
            log_warn(paste(
                "lipidomics artifact read-through failed without changing state:",
                conditionMessage(error)
            ))
            list(
                enabled = artifact_project,
                ok = FALSE,
                resumed = FALSE,
                reason = "artifact_readthrough_failed",
                error_class = class(error),
                error_message = conditionMessage(error)
            )
        }
    )
    result$artifact_project <- artifact_project
    recordArtifactStageResult(workflow_data, "resume", result)
}

previewLipidImportArtifact <- function(
    workflow_data,
    assay_name,
    projections = NULL,
    limit = 100L
) {
    if (!workflowCapabilityScalarString(assay_name) ||
        length(limit) != 1L || !is.numeric(limit) || is.na(limit) ||
        limit < 1L || limit > 1000L || limit != floor(limit)) {
        lipidArtifactAbort(
            "lipidomics preview requires one assay and a limit from 1 to 1000",
            "multischolar_invalid_lipidomics_preview"
        )
    }
    context <- workflow_data$workflow_context
    refs <- workflow_data$artifact_readthrough_refs$import
    proof <- workflow_data$artifact_readthrough_proof
    valid <- inherits(context, "WorkflowContext") && context$isBound() &&
        is.list(refs) && is.list(proof) &&
        identical(proof$descriptor_id,
            artifactLipidomicsWorkflowDescriptor()$descriptor_id)
    if (!isTRUE(valid)) {
        lipidArtifactAbort(
            "lipidomics preview requires a resumed coordinator",
            "multischolar_invalid_lipidomics_preview"
        )
    }
    assay_order <- proof$assay_order
    roles <- sprintf("assay_%04d", seq_along(assay_order))
    index <- match(assay_name, assay_order)
    if (is.na(index)) {
        lipidArtifactAbort(
            "lipidomics preview assay is not in the committed import",
            "multischolar_invalid_lipidomics_preview"
        )
    }
    store <- newArtifactStore(
        context$getPaths(),
        context$getIdentity()$project_id
    )
    value <- artifactSettledResumeReadLocator(
        lipidArtifactReadthroughAdapter(),
        store,
        refs[[roles[[index]]]]
    )
    if (!is.null(projections)) {
        if (!is.character(projections) || any(!projections %in% names(value))) {
            lipidArtifactAbort(
                "lipidomics preview projection is not in the assay schema",
                "multischolar_invalid_lipidomics_preview"
            )
        }
        value <- value[projections]
    }
    utils::head(value, as.integer(limit))
}
