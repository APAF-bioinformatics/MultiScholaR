# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

.METAB_READTHROUGH_PROOF_SCHEMA <- "multischolar.metabolomics.readthrough_proof"
.METAB_READTHROUGH_PROOF_VERSION <- 1L

metabArtifactContractVersion <- function(value) {
    artifactStageReadthroughContractVersion(
        value,
        .METAB_READTHROUGH_VERSION
    )
}

validateMetabResumeContract <- function(evidence) {
    adapter <- metabArtifactReadthroughAdapter()
    import <- evidence$import
    design <- evidence$design
    import_parameters <- import$parameters
    design_parameters <- design$parameters
    import_order <- unlist(import_parameters$assay_order, use.names = FALSE)
    design_order <- unlist(design_parameters$assay_order, use.names = FALSE)
    valid <- metabArtifactContractVersion(
        import_parameters$readthrough_contract_version
    ) && metabArtifactContractVersion(
        design_parameters$readthrough_contract_version
    ) && identical(
        import_parameters$capability_id,
        evidence$descriptor$descriptor_id
    ) && identical(
        design_parameters$capability_id,
        evidence$descriptor$descriptor_id
    ) && identical(import_parameters$input_format, "custom") &&
        identical(import_parameters$data_level, "metabolite") &&
        identical(design_parameters$state_name, "metab_raw_data_s4") &&
        identical(design_parameters$workflow_type, "metabolomics_standard") &&
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
            "metabolomics import, design, and state are not one exact lineage"
        )
    }
    invisible(TRUE)
}

newMetabResumeStateManager <- function(context, resource_policy = NULL) {
    ArtifactWorkflowState$new(
        workflow_context = context,
        resource_policy = resource_policy,
        dehydrate_fn = dehydrateMetabolomicsS4Artifact,
        validate_bundle_fn = validateMetabolomicsS4Bundle,
        hydrate_fn = hydrateMetabolomicsS4Artifact,
        descriptor_contract = artifactStageDescriptorContract(
            artifactMetabolomicsWorkflowDescriptor()
        )
    )
}

validateMetabStateEvidence <- function(manager, evidence, state_object) {
    adapter <- metabArtifactReadthroughAdapter()
    metadata <- manager$getStateMetadata()
    audit <- metadata$audit_metadata
    expected_refs <- artifactStageNormalizedRefs(evidence$design$refs)
    actual_refs <- artifactStageNormalizedRefs(audit$stage_artifact_refs)
    valid <- identical(
        manager$getCurrentGenerationId(),
        evidence$current_state$generation_id
    ) && identical(
        manager$getCurrentStateName(),
        evidence$design$parameters$state_name
    ) && identical(audit$stage_id, "design") &&
        identical(audit$run_id, evidence$design$run_id) &&
        metabArtifactContractVersion(audit$readthrough_contract_version) &&
        identical(audit$parent_import_run_id, evidence$import$run_id) &&
        identical(
            audit$parent_import_generation_id,
            evidence$import$generation_id
        ) && identical(actual_refs, expected_refs) &&
        methods::is(state_object, "MetaboliteAssayData") &&
        identical(methods::validObject(state_object, test = TRUE), TRUE)
    if (!isTRUE(valid)) {
        artifactStageReadthroughAbort(
            adapter,
            "incomplete_contract",
            "metabolomics current S4 state differs from its design generation"
        )
    }
    invisible(TRUE)
}

metabResumeAssayTables <- function(tables, parameters, prefix = "assay") {
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

validateMetabImportTables <- function(evidence, tables, mapping) {
    adapter <- metabArtifactReadthroughAdapter()
    parameters <- evidence$import$parameters
    assays <- metabResumeAssayTables(tables, parameters)
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
            "metabolomics import tables differ from their assay manifest"
        )
    }
    assays
}

validateMetabDesignTables <- function(
    evidence,
    tables,
    mapping,
    state_object,
    config
) {
    adapter <- metabArtifactReadthroughAdapter()
    parameters <- evidence$design$parameters
    assays <- metabResumeAssayTables(tables, parameters, "cleaned")
    alignment <- tables$assay_alignment
    dependencies <- tables$raw_s4_dependencies
    state_assays <- methods::slot(state_object, "metabolite_data")
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
            feature_column <- methods::slot(state_object, "metabolite_id_column")
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
            "metabolomics design tables differ from the current S4 state"
        )
    }
    assays
}

hydrateMetabResumeBundle <- function(
    context,
    resource_policy = NULL,
    retain_source_payloads = TRUE
) {
    adapter <- metabArtifactReadthroughAdapter()
    if (!identical(retain_source_payloads, TRUE)) {
        artifactStageReadthroughAbort(
            adapter,
            "incomplete_contract",
            "metabolomics import/design payloads must remain retained"
        )
    }
    evidence <- collectMetabResumeEvidence(
        context,
        resource_policy,
        payload_validation = "sidecar"
    )
    validateMetabResumeContract(evidence)
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
    manager <- newMetabResumeStateManager(context, resource_policy)
    manager_owned <- TRUE
    on.exit({
        if (manager_owned) try(manager$close(), silent = TRUE)
    }, add = TRUE)
    state_object <- manager$getState()
    validateMetabStateEvidence(manager, evidence, state_object)
    config <- manager$getStateConfig()
    mapping <- artifactStageColumnMapping(
        adapter,
        evidence$import$parameters
    )
    data_tbl <- validateMetabImportTables(evidence, import_tables, mapping)
    data_cln <- validateMetabDesignTables(
        evidence,
        design_tables,
        mapping,
        state_object,
        config
    )
    contrasts <- artifactStageRestoreContrasts(
        adapter,
        design_tables$contrasts,
        evidence$design$parameters$contrasts_kind
    )
    preflight <- validateMetabDesignDaPreflight(
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
            "metabolomics resumed design fails scientific preflight"
        )
    }
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
        source_payloads_retained = TRUE,
        readthrough_mode = "full_retained",
        evidence = evidence
    )
    manager_owned <- FALSE
    bundle
}

metabReadthroughConsumerInventory <- function() {
    readers <- c(
        "previewMetabImportArtifact", "validateMetabDesignDaPreflight",
        "runMetabQcS4ServerBody", "normaliseBetweenSamples",
        "buildMetabNormExportSessionData", "differentialAbundanceAnalysisHelper",
        "getMetabSummaryStateObject"
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
            envir = environment(metabReadthroughConsumerInventory),
            mode = "function",
            inherits = TRUE
        ), logical(1)),
        stringsAsFactors = FALSE
    )
}

validateMetabReadthroughConsumerInventory <- function(inventory) {
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
        metabArtifactAbort(
            "metabolomics read-through consumer inventory is incomplete",
            "multischolar_incomplete_metabolomics_consumer_inventory"
        )
    }
    inventory
}

metabReadthroughResourceSnapshot <- function(workflow_data, phase) {
    fields <- c("data_tbl", "data_cln", "design_matrix")
    object_bytes <- vapply(fields, \(field) {
        value <- tryCatch(workflow_data[[field]], error = \(...) NULL)
        unname(as.numeric(utils::object.size(value)))
    }, numeric(1))
    manager <- tryCatch(workflow_data$state_manager, error = \(...) NULL)
    state <- tryCatch(manager$getState(), error = \(...) NULL)
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

newMetabReadthroughProof <- function(bundle, inventory) {
    evidence <- bundle$evidence
    list(
        schema = .METAB_READTHROUGH_PROOF_SCHEMA,
        schema_version = .METAB_READTHROUGH_PROOF_VERSION,
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
        source_payloads_retained = TRUE,
        eviction_performed = FALSE,
        consumer_inventory_digest = artifactSemanticDigest(inventory),
        verified_at = artifactRefUtcNow()
    )
}

captureMetabResumeFields <- function(workflow_data, fields) {
    values <- lapply(fields, \(field) {
        tryCatch(workflow_data[[field]], error = \(...) NULL)
    })
    stats::setNames(values, fields)
}

restoreMetabResumeFields <- function(workflow_data, values) {
    for (field in names(values)) workflow_data[[field]] <- values[[field]]
    invisible(TRUE)
}

applyMetabResumeBundle <- function(
    workflow_data,
    bundle,
    failure_injector = NULL,
    inventory_fn = metabReadthroughConsumerInventory
) {
    inventory <- validateMetabReadthroughConsumerInventory(inventory_fn())
    fields <- c(
        "workflow_context", "state_manager", "data_tbl", "data_cln",
        "design_matrix", "contrasts_tbl", "config_list", "column_mapping",
        "data_format", "data_type", "artifact_readthrough_refs",
        "artifact_readthrough_proof", "artifact_compatibility_checkpoint",
        "tab_status"
    )
    previous <- captureMetabResumeFields(workflow_data, fields)
    committed <- FALSE
    on.exit({
        if (!committed) restoreMetabResumeFields(workflow_data, previous)
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
    workflow_data$data_format <- "custom"
    workflow_data$data_type <- "metabolite"
    workflow_data$artifact_readthrough_refs <- list(
        import = bundle$evidence$import$refs,
        design = bundle$evidence$design$refs
    )
    proof <- newMetabReadthroughProof(bundle, inventory)
    checkpoint <- list(
        strategy = "restore_memory_reads_keep_metabolomics_generations",
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

resumeMetabArtifactWorkflow <- function(
    workflow_data,
    experiment_paths,
    experiment_label,
    storage_policy = NULL,
    resource_policy = NULL,
    failure_injector = NULL,
    inventory_fn = metabReadthroughConsumerInventory
) {
    prepared <- createMetabResumeContext(
        experiment_paths,
        experiment_label,
        storage_policy
    )
    if (!isTRUE(prepared$enabled)) {
        return(c(prepared, list(ok = TRUE, resumed = FALSE)))
    }
    before <- metabReadthroughResourceSnapshot(workflow_data, "before")
    bundle <- hydrateMetabResumeBundle(
        prepared$context,
        resource_policy,
        retain_source_payloads = TRUE
    )
    applied <- FALSE
    on.exit({
        if (!applied) try(bundle$state_manager$close(), silent = TRUE)
    }, add = TRUE)
    readthrough <- applyMetabResumeBundle(
        workflow_data,
        bundle,
        failure_injector,
        inventory_fn
    )
    applied <- TRUE
    after <- metabReadthroughResourceSnapshot(workflow_data, "after")
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
        assay_order = names(bundle$data_cln),
        source_payloads_retained = TRUE,
        eviction_performed = FALSE,
        compatibility_checkpoint = readthrough$compatibility_checkpoint,
        resource_evidence = list(before = before, after = after)
    )
}

resumeMetabArtifactWorkflowSafely <- function(
    workflow_data,
    experiment_paths,
    experiment_label,
    storage_policy = NULL,
    resource_policy = NULL,
    failure_injector = NULL,
    inventory_fn = metabReadthroughConsumerInventory,
    log_warn = logger::log_warn
) {
    artifact_project <- metabArtifactProjectPresent(
        experiment_paths,
        experiment_label,
        storage_policy
    )
    result <- tryCatch(
        resumeMetabArtifactWorkflow(
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
                "metabolomics artifact read-through failed without changing state:",
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

previewMetabImportArtifact <- function(
    workflow_data,
    assay_name,
    projections = NULL,
    limit = 100L
) {
    if (!workflowCapabilityScalarString(assay_name) ||
        length(limit) != 1L || !is.numeric(limit) || is.na(limit) ||
        limit < 1L || limit > 1000L || limit != floor(limit)) {
        metabArtifactAbort(
            "metabolomics preview requires one assay and a limit from 1 to 1000",
            "multischolar_invalid_metabolomics_preview"
        )
    }
    context <- workflow_data$workflow_context
    refs <- workflow_data$artifact_readthrough_refs$import
    proof <- workflow_data$artifact_readthrough_proof
    valid <- inherits(context, "WorkflowContext") && context$isBound() &&
        is.list(refs) && is.list(proof) &&
        identical(proof$descriptor_id,
            artifactMetabolomicsWorkflowDescriptor()$descriptor_id)
    if (!isTRUE(valid)) {
        metabArtifactAbort(
            "metabolomics preview requires a resumed coordinator",
            "multischolar_invalid_metabolomics_preview"
        )
    }
    assay_order <- proof$assay_order
    roles <- sprintf("assay_%04d", seq_along(assay_order))
    index <- match(assay_name, assay_order)
    if (is.na(index)) {
        metabArtifactAbort(
            "metabolomics preview assay is not in the committed import",
            "multischolar_invalid_metabolomics_preview"
        )
    }
    store <- newArtifactStore(
        context$getPaths(),
        context$getIdentity()$project_id
    )
    value <- artifactStageReadTable(
        metabArtifactReadthroughAdapter(),
        store,
        refs[[roles[[index]]]]
    )
    if (!is.null(projections)) {
        if (!is.character(projections) || any(!projections %in% names(value))) {
            metabArtifactAbort(
                "metabolomics preview projection is not in the assay schema",
                "multischolar_invalid_metabolomics_preview"
            )
        }
        value <- value[projections]
    }
    utils::head(value, as.integer(limit))
}
