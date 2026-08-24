# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

.PROT_NONDIA_READTHROUGH_PROOF_SCHEMA <-
    "multischolar.prot_nondia.readthrough_proof"
.PROT_NONDIA_READTHROUGH_PROOF_VERSION <- 1L

#' Test a non-DIA read-through contract version
#' @param value Serialized version value.
#' @return A scalar logical.
#' @noRd
protNonDiaArtifactContractVersion <- function(value) {
    artifactStageReadthroughContractVersion(
        value,
        .PROT_NONDIA_READTHROUGH_VERSION
    )
}

#' Validate parent-linked non-DIA import, design, and state evidence
#' @param evidence Collected read-through evidence.
#' @return `TRUE`, invisibly.
#' @noRd
protNonDiaValidateResumeContract <- function(evidence) {
    descriptor <- protNonDiaReadthroughDescriptor(
        evidence$descriptor$descriptor_id
    )
    adapter <- protNonDiaArtifactReadthroughAdapter(descriptor)
    spec <- protNonDiaArtifactImportSpec(descriptor$identity$input_format)
    design <- evidence$design
    import <- evidence$import
    design_ref <- design$refs$cleaned_data
    import_ref <- import$refs$canonical_data
    design_parameters <- design$parameters
    import_parameters <- import$parameters
    mapping <- artifactStageColumnMapping(adapter, import_parameters)
    valid <- protNonDiaArtifactContractVersion(
        design_parameters$readthrough_contract_version
    ) && protNonDiaArtifactContractVersion(
        import_parameters$readthrough_contract_version
    ) && identical(design_parameters$capability_id, descriptor$descriptor_id) &&
        identical(import_parameters$capability_id, descriptor$descriptor_id) &&
        identical(design_parameters$workflow_type, spec$workflow_type) &&
        identical(import_parameters$workflow_type, spec$workflow_type) &&
        identical(import_parameters$input_format, descriptor$identity$input_format) &&
        identical(import_parameters$data_level, "protein") &&
        identical(import_parameters$parser_parameters, spec$parser_parameters) &&
        identical(import_parameters$column_mapping, mapping) &&
        identical(design_parameters$state_name, "protein_s4_initial") &&
        identical(design_parameters$parent_import_run_id, import$run_id) &&
        identical(
            design_parameters$parent_import_generation_id,
            import$generation_id
        ) && identical(
            design_parameters$parent_import_artifact_id,
            import_ref$artifact_id
        ) && identical(
            design_parameters$parent_import_semantic_digest,
            import_ref$hash_policy$semantic$digest
        ) && identical(
            evidence$current_state$logical_name,
            design_parameters$state_name
        ) && identical(
            design_ref$logical_key$generation_id,
            design$generation_id
        ) && artifactStageIdentityMatches(evidence$identity, descriptor)
    if (!isTRUE(valid)) {
        artifactStageReadthroughAbort(
            adapter,
            "incomplete_contract",
            "non-DIA import, design, and state evidence is not one lineage",
            project_id = evidence$identity$project_id,
            design_run_id = design$run_id,
            import_run_id = import$run_id
        )
    }
    invisible(TRUE)
}

#' Validate the current non-DIA protein S4 state and lineage
#' @param manager Descriptor-bound artifact state manager.
#' @param evidence Collected read-through evidence.
#' @param state_object Independently hydrated S4 state.
#' @return `TRUE`, invisibly.
#' @noRd
protNonDiaValidateStateEvidence <- function(manager, evidence, state_object) {
    descriptor <- evidence$descriptor
    adapter <- protNonDiaArtifactReadthroughAdapter(descriptor)
    metadata <- manager$getStateMetadata()
    audit <- metadata$audit_metadata
    expected_refs <- artifactStageNormalizedRefs(evidence$design$refs)
    actual_refs <- artifactStageNormalizedRefs(audit$stage_artifact_refs)
    role <- artifactProteomicsNonDiaCodecRole(
        descriptor$descriptor_id,
        evidence$design$parameters$state_name,
        descriptor$descriptor_version
    )
    shape_valid <- tryCatch(
        {
            artifactValidateProteomicsNonDiaProteinState(
                state_object,
                role,
                evidence$design$parameters$state_name
            )
            TRUE
        },
        error = \(...) FALSE
    )
    valid <- identical(
        manager$getCurrentGenerationId(),
        evidence$current_state$generation_id
    ) && identical(
        manager$getCurrentStateName(),
        evidence$design$parameters$state_name
    ) && identical(audit$stage_id, "design") &&
        identical(audit$run_id, evidence$design$run_id) &&
        protNonDiaArtifactContractVersion(audit$readthrough_contract_version) &&
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
        ) && identical(actual_refs, expected_refs) && isTRUE(shape_valid) &&
        isS4(state_object) &&
        identical(class(state_object)[[1L]], "ProteinQuantitativeData") &&
        identical(methods::validObject(state_object, test = TRUE), TRUE)
    if (!isTRUE(valid)) {
        artifactStageReadthroughAbort(
            adapter,
            "incomplete_contract",
            "non-DIA current S4 state does not match its committed design run",
            design_run_id = evidence$design$run_id,
            state_generation_id = evidence$current_state$generation_id
        )
    }
    invisible(TRUE)
}

#' Validate non-DIA hydrated scientific tables against the S4 state
#' @param evidence Collected read-through evidence.
#' @param state_object Independently hydrated S4 state.
#' @param config Independently hydrated state configuration.
#' @param import_tables Hydrated import tables.
#' @param design_tables Hydrated design tables.
#' @return `TRUE`, invisibly.
#' @noRd
protNonDiaValidateScientificTables <- function(
    evidence,
    state_object,
    config,
    import_tables,
    design_tables
) {
    adapter <- protNonDiaArtifactReadthroughAdapter(evidence$descriptor)
    expected_args <- artifactStageMetadataTable(state_object@args)
    valid <- identical(
        import_tables$canonical_data,
        design_tables$cleaned_data
    ) && identical(
        design_tables$design_matrix,
        state_object@design_matrix
    ) && identical(design_tables$args, expected_args) &&
        identical(config, state_object@args)
    if (!isTRUE(valid)) {
        artifactStageReadthroughAbort(
            adapter,
            "incomplete_contract",
            "non-DIA design tables differ from the current scientific S4 state"
        )
    }
    invisible(TRUE)
}

#' Validate settled non-DIA metadata against the current S4 state
#' @param evidence Collected read-through evidence.
#' @param state_object Independently hydrated S4 state.
#' @param config Independently hydrated state configuration.
#' @param design_tables Hydrated small design tables.
#' @return `TRUE`, invisibly.
#' @noRd
protNonDiaValidateSettledScientificTables <- function(
    evidence,
    state_object,
    config,
    design_tables
) {
    adapter <- protNonDiaArtifactReadthroughAdapter(evidence$descriptor)
    expected_args <- artifactStageMetadataTable(state_object@args)
    valid <- identical(
        design_tables$design_matrix,
        state_object@design_matrix
    ) && identical(design_tables$args, expected_args) &&
        identical(config, state_object@args)
    if (!isTRUE(valid)) {
        artifactStageReadthroughAbort(
            adapter,
            "incomplete_contract",
            "non-DIA settled metadata differs from the current scientific S4 state"
        )
    }
    invisible(TRUE)
}

#' Construct a descriptor-bound non-DIA artifact state manager
#' @param context Exact bound workflow context.
#' @param descriptor Exact supported workflow descriptor.
#' @param resource_policy Optional project registry resource policy.
#' @return An artifact-backed workflow state manager.
#' @noRd
newProtNonDiaResumeStateManager <- function(
    context,
    descriptor,
    resource_policy = NULL
) {
    newWorkflowState(
        workflow_context = context,
        resource_policy = resource_policy,
        workflow_descriptor = descriptor,
        descriptor_catalogue = artifactWorkflowDescriptorCatalogue(),
        codec_catalogue = artifactS4CodecCatalogue()
    )
}

#' Hydrate a complete non-DIA proteomics resume bundle
#' @param context Exact bound workflow context.
#' @param descriptor Exact supported workflow descriptor.
#' @param resource_policy Optional project registry resource policy.
#' @param retain_source_payloads Must remain `TRUE` until the eviction ticket.
#' @return A fully validated resume bundle that owns an open state manager.
#' @noRd
hydrateProtNonDiaResumeBundle <- function(
    context,
    descriptor,
    resource_policy = NULL,
    retain_source_payloads = TRUE
) {
    adapter <- protNonDiaArtifactReadthroughAdapter(descriptor)
    if (!is.logical(retain_source_payloads) ||
        length(retain_source_payloads) != 1L || is.na(retain_source_payloads)) {
        artifactStageReadthroughAbort(
            adapter,
            "incomplete_contract",
            "non-DIA source retention mode must be true or false"
        )
    }
    evidence <- collectProtNonDiaResumeEvidence(
        context,
        descriptor,
        resource_policy,
        payload_validation = if (isTRUE(retain_source_payloads)) {
            "sidecar"
        } else {
            "digest"
        }
    )
    protNonDiaValidateResumeContract(evidence)
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
    } else {
        import_tables <- list(canonical_data = NULL)
        design_tables <- readArtifactStageRoles(
            adapter,
            evidence$store,
            evidence$design,
            c("design_matrix", "contrasts", "args", "annotations", "sequences")
        )
        design_tables$cleaned_data <- NULL
    }
    manager <- newProtNonDiaResumeStateManager(
        context,
        descriptor,
        resource_policy
    )
    manager_owned <- TRUE
    on.exit({
        if (manager_owned) try(manager$close(), silent = TRUE)
    }, add = TRUE)
    state_object <- manager$getState()
    protNonDiaValidateStateEvidence(manager, evidence, state_object)
    config <- manager$getStateConfig()
    if (isTRUE(retain_source_payloads)) {
        protNonDiaValidateScientificTables(
            evidence,
            state_object,
            config,
            import_tables,
            design_tables
        )
    } else {
        protNonDiaValidateSettledScientificTables(
            evidence,
            state_object,
            config,
            design_tables
        )
    }
    parameters <- evidence$design$parameters
    bundle <- list(
        context = context,
        descriptor = descriptor,
        state_manager = manager,
        state_object = state_object,
        data_tbl = import_tables$canonical_data,
        data_cln = design_tables$cleaned_data,
        design_matrix = design_tables$design_matrix,
        contrasts_tbl = artifactStageRestoreContrasts(
            adapter,
            design_tables$contrasts,
            parameters$contrasts_kind
        ),
        config_list = config,
        column_mapping = artifactStageColumnMapping(
            adapter,
            evidence$import$parameters
        ),
        uniprot_dat_cln = artifactStageRestoreOptional(
            adapter,
            design_tables$annotations,
            parameters$annotation_available,
            "annotations"
        ),
        aa_seq_tbl_final = artifactStageRestoreOptional(
            adapter,
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

#' Inventory every non-DIA import/design read-through consumer class
#' @return A data frame of consumer readers and their restored source.
#' @noRd
protNonDiaReadthroughConsumerInventory <- function() {
    readers <- c(
        "previewProtNonDiaImportArtifact",
        "hydrateArtifactStageInput",
        "normaliseBetweenSamples",
        "collectProtNormExportSessionData",
        "resolveProtDaAvailableContrasts",
        "runTestsContrasts",
        "resolveProtSummaryExpectedTemplate",
        "prepareProtSummaryCopyInputs"
    )
    data.frame(
        category = c(
            "preview", "normalization", "normalization", "session", "da",
            "da", "report", "compatibility"
        ),
        reader = readers,
        restored_source = c(
            "bounded immutable import ref",
            "current artifact state generation",
            "hydrated protein S4 state",
            "current state plus restored coordinator metadata",
            "restored contrasts compatibility adapter",
            "hydrated quantitative and design matrices",
            "restored LFQ or TMT workflow type",
            "restored design and contrasts tables"
        ),
        verified = vapply(readers, \(reader) {
            exists(
                reader,
                envir = environment(protNonDiaReadthroughConsumerInventory),
                mode = "function",
                inherits = TRUE
            )
        }, logical(1)),
        stringsAsFactors = FALSE
    )
}

#' Validate the non-DIA read-through consumer inventory
#' @param inventory Candidate consumer inventory.
#' @return The validated inventory.
#' @noRd
validateProtNonDiaReadthroughConsumerInventory <- function(inventory) {
    required_categories <- c(
        "preview", "normalization", "session", "da", "report",
        "compatibility"
    )
    valid <- is.data.frame(inventory) && identical(
        names(inventory),
        c("category", "reader", "restored_source", "verified")
    ) && setequal(unique(inventory$category), required_categories) &&
        !anyDuplicated(inventory$reader) && all(inventory$verified) &&
        all(nzchar(inventory$reader)) && all(nzchar(inventory$restored_source))
    if (!isTRUE(valid)) {
        protNonDiaArtifactAbort(
            "non-DIA read-through consumer inventory is incomplete",
            "multischolar_incomplete_prot_nondia_consumer_inventory"
        )
    }
    inventory
}

#' Build a payload-free proof for non-DIA read-through
#' @param bundle Fully validated non-DIA resume bundle.
#' @return A payload-free read-through proof.
#' @noRd
newProtNonDiaReadthroughProof <- function(bundle) {
    evidence <- bundle$evidence
    list(
        schema = .PROT_NONDIA_READTHROUGH_PROOF_SCHEMA,
        schema_version = .PROT_NONDIA_READTHROUGH_PROOF_VERSION,
        descriptor_id = evidence$descriptor$descriptor_id,
        descriptor_version = evidence$descriptor$descriptor_version,
        project_id = evidence$identity$project_id,
        workflow_id = evidence$identity$workflow_id,
        import_run_id = evidence$import$run_id,
        design_run_id = evidence$design$run_id,
        state_generation_id = evidence$current_state$generation_id,
        payload_validated = TRUE,
        registry_ready = TRUE,
        current_pointer_verified = TRUE,
        readthrough_completed = TRUE,
        annotation_completed = isTRUE(bundle$annotation_completed),
        readthrough_mode = bundle$readthrough_mode,
        source_payloads_retained = bundle$source_payloads_retained,
        import_refs = artifactPayloadReadthroughRefProof(evidence$import$refs),
        design_refs = artifactPayloadReadthroughRefProof(evidence$design$refs),
        consumer_inventory_digest = artifactSemanticDigest(
            protNonDiaReadthroughConsumerInventory()
        ),
        verified_at = artifactRefUtcNow()
    )
}

#' Build the rollback checkpoint for non-DIA read-through
#' @param proof Validated non-DIA read-through proof.
#' @return A payload-free compatibility checkpoint.
#' @noRd
newProtNonDiaCompatibilityCheckpoint <- function(proof) {
    list(
        strategy = "restore_memory_reads_keep_artifact_generations",
        verified = TRUE,
        descriptor_id = proof$descriptor_id,
        project_id = proof$project_id,
        workflow_id = proof$workflow_id,
        design_run_id = proof$design_run_id,
        state_generation_id = proof$state_generation_id,
        consumer_inventory_digest = proof$consumer_inventory_digest,
        verified_at = proof$verified_at
    )
}

#' Capture tuple-specific read-through memory characterization
#' @param workflow_data Mutable proteomics workflow state.
#' @param phase Characterization phase label.
#' @return RSS and retained workflow object sizes without invoking GC.
#' @noRd
protNonDiaReadthroughResourceSnapshot <- function(workflow_data, phase) {
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

#' Capture workflow fields used by atomic non-DIA resume application
#' @param workflow_data Mutable proteomics workflow state.
#' @param fields Exact fields to capture.
#' @return A named list of prior field values.
#' @noRd
captureProtNonDiaResumeFields <- function(workflow_data, fields) {
    values <- lapply(fields, \(field) {
        tryCatch(workflow_data[[field]], error = \(...) NULL)
    })
    stats::setNames(values, fields)
}

#' Restore workflow fields after a failed non-DIA resume apply
#' @param workflow_data Mutable proteomics workflow state.
#' @param values Named prior field values.
#' @return `TRUE`, invisibly.
#' @noRd
restoreProtNonDiaResumeFields <- function(workflow_data, values) {
    for (field in names(values)) workflow_data[[field]] <- values[[field]]
    invisible(TRUE)
}

#' Apply a validated non-DIA resume bundle atomically
#' @param workflow_data Mutable proteomics workflow state.
#' @param bundle Fully validated non-DIA resume bundle.
#' @param failure_injector Optional failure injector used by tests.
#' @param inventory_fn Consumer inventory provider.
#' @return Read-through proof and compatibility checkpoint, invisibly.
#' @noRd
applyProtNonDiaResumeBundle <- function(
    workflow_data,
    bundle,
    failure_injector = NULL,
    inventory_fn = protNonDiaReadthroughConsumerInventory
) {
    inventory <- validateProtNonDiaReadthroughConsumerInventory(inventory_fn())
    fields <- c(
        "workflow_context", "state_manager", "data_tbl", "data_cln",
        "design_matrix", "contrasts_tbl", "config_list", "column_mapping",
        "uniprot_dat_cln", "aa_seq_tbl_final", "data_format", "data_type",
        "artifact_readthrough_refs", "artifact_readthrough_proof",
        "artifact_compatibility_checkpoint", "tab_status"
    )
    previous <- captureProtNonDiaResumeFields(workflow_data, fields)
    dependencies <- c(
        "design_matrix", "contrasts", "config", "annotations", "sequences"
    )
    previous_globals <- captureProtContextLegacyGlobals(dependencies)
    committed <- FALSE
    on.exit({
        if (!committed) {
            restoreProtNonDiaResumeFields(workflow_data, previous)
            restoreProtContextLegacyGlobals(previous_globals)
        }
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
    workflow_data$uniprot_dat_cln <- bundle$uniprot_dat_cln
    workflow_data$aa_seq_tbl_final <- bundle$aa_seq_tbl_final
    workflow_data$data_format <- bundle$descriptor$identity$input_format
    workflow_data$data_type <- "protein"
    workflow_data$artifact_readthrough_refs <- list(
        import = bundle$evidence$import$refs,
        design = bundle$evidence$design$refs
    )
    proof <- newProtNonDiaReadthroughProof(bundle)
    proof$consumer_inventory_digest <- artifactSemanticDigest(inventory)
    checkpoint <- newProtNonDiaCompatibilityCheckpoint(proof)
    workflow_data$artifact_readthrough_proof <- proof
    workflow_data$artifact_compatibility_checkpoint <- checkpoint
    status <- workflow_data$tab_status
    status$setup_import <- "complete"
    status$design_matrix <- "complete"
    status$quality_control <- "complete"
    if (identical(status$normalization, "disabled")) {
        status$normalization <- "pending"
    }
    workflow_data$tab_status <- status
    publishProtContextLegacyGlobal(
        "design_matrix", bundle$design_matrix, workflow_data = workflow_data
    )
    publishProtContextLegacyGlobal(
        "contrasts", bundle$contrasts_tbl, workflow_data = workflow_data
    )
    publishProtContextLegacyGlobal(
        "config", bundle$config_list, workflow_data = workflow_data
    )
    publishProtContextLegacyGlobal(
        "annotations", bundle$uniprot_dat_cln, workflow_data = workflow_data
    )
    publishProtContextLegacyGlobal(
        "sequences", bundle$aa_seq_tbl_final, workflow_data = workflow_data
    )
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

#' Test whether a non-DIA resume error can use an explicit legacy adapter
#' @param error Candidate resume error.
#' @return A scalar logical.
#' @noRd
protNonDiaArtifactLegacyEligibleError <- function(error) {
    inherits(error, c(
        "multischolar_missing_prot_nondia_committed_stage",
        "multischolar_incomplete_prot_nondia_readthrough_contract",
        "multischolar_missing_prot_nondia_current_state"
    ))
}

#' Resume one exact non-DIA proteomics artifact workflow
#' @param workflow_data Mutable proteomics workflow state.
#' @param experiment_paths Workflow project paths.
#' @param experiment_label Project experiment label.
#' @param capability_id Optional exact workflow capability identifier.
#' @param storage_policy Optional workflow storage policy.
#' @param resource_policy Optional project registry resource policy.
#' @param legacy_adapter Optional whole-workflow memory migration adapter.
#' @param failure_injector Optional failure injector used by tests.
#' @param inventory_fn Consumer inventory provider.
#' @return A structured resume result.
#' @noRd
resumeProtNonDiaArtifactWorkflow <- function(
    workflow_data,
    experiment_paths,
    experiment_label,
    capability_id = NULL,
    storage_policy = NULL,
    resource_policy = NULL,
    legacy_adapter = NULL,
    failure_injector = NULL,
    inventory_fn = protNonDiaReadthroughConsumerInventory
) {
    prepared <- createProtNonDiaResumeContext(
        experiment_paths,
        experiment_label,
        capability_id,
        storage_policy
    )
    if (!isTRUE(prepared$enabled)) {
        return(c(prepared, list(ok = TRUE, resumed = FALSE)))
    }
    descriptor <- prepared$descriptor
    before <- protNonDiaReadthroughResourceSnapshot(workflow_data, "before")
    retain_source_payloads <- !identical(
        prepared$context$getStorageDecision()$effective_rollout,
        "evict"
    )
    bundle <- tryCatch(
        hydrateProtNonDiaResumeBundle(
            prepared$context,
            descriptor,
            resource_policy,
            retain_source_payloads = retain_source_payloads
        ),
        error = \(error) error
    )
    if (inherits(bundle, "error")) {
        if (protNonDiaArtifactLegacyEligibleError(bundle) &&
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
    readthrough <- applyProtNonDiaResumeBundle(
        workflow_data,
        bundle,
        failure_injector,
        inventory_fn
    )
    applied <- TRUE
    after <- protNonDiaReadthroughResourceSnapshot(workflow_data, "after")
    list(
        enabled = TRUE,
        ok = TRUE,
        resumed = TRUE,
        reason = "artifact_readthrough",
        capability_id = descriptor$descriptor_id,
        project_id = bundle$evidence$identity$project_id,
        import_run_id = bundle$evidence$import$run_id,
        design_run_id = bundle$evidence$design$run_id,
        state_generation_id = bundle$evidence$current_state$generation_id,
        import_ref = unclass(bundle$evidence$import$refs$canonical_data),
        source_payloads_retained = bundle$source_payloads_retained,
        readthrough_mode = bundle$readthrough_mode,
        compatibility_checkpoint = readthrough$compatibility_checkpoint,
        resource_evidence = list(before = before, after = after)
    )
}

#' Safely resume one exact non-DIA proteomics artifact workflow
#' @param workflow_data Mutable proteomics workflow state.
#' @param experiment_paths Workflow project paths.
#' @param experiment_label Project experiment label.
#' @param capability_id Optional exact workflow capability identifier.
#' @param storage_policy Optional workflow storage policy.
#' @param resource_policy Optional project registry resource policy.
#' @param failure_injector Optional failure injector used by tests.
#' @param inventory_fn Consumer inventory provider.
#' @param log_warn Warning logger used for additive resume failures.
#' @return The recorded resume result, invisibly.
#' @noRd
resumeProtNonDiaArtifactWorkflowSafely <- function(
    workflow_data,
    experiment_paths,
    experiment_label,
    capability_id = NULL,
    storage_policy = NULL,
    resource_policy = NULL,
    failure_injector = NULL,
    inventory_fn = protNonDiaReadthroughConsumerInventory,
    log_warn = logger::log_warn
) {
    artifact_project <- protNonDiaArtifactProjectPresent(
        experiment_paths,
        experiment_label,
        storage_policy
    )
    result <- tryCatch(
        resumeProtNonDiaArtifactWorkflow(
            workflow_data,
            experiment_paths,
            experiment_label,
            capability_id,
            storage_policy,
            resource_policy,
            failure_injector = failure_injector,
            inventory_fn = inventory_fn
        ),
        error = \(error) {
            log_warn(paste(
                "non-DIA artifact read-through failed without changing workflow state:",
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
    recordArtifactStageResult(workflow_data, "resume", result)
}

#' Query a bounded non-DIA import preview
#' @param context Exact bound workflow context.
#' @param import_ref Optional immutable import artifact reference.
#' @param projections Optional allowlisted projections.
#' @param filters Optional bounded filters.
#' @param limit Optional bounded row limit.
#' @param resource_policy Optional project registry resource policy.
#' @param query_session Optional reusable bounded query session.
#' @return A bounded import preview data frame.
#' @noRd
previewProtNonDiaImportArtifact <- function(
    context,
    import_ref = NULL,
    projections = NULL,
    filters = list(),
    limit = NULL,
    resource_policy = NULL,
    query_session = NULL
) {
    if (!inherits(context, "WorkflowContext") || !context$isBound()) {
        protNonDiaArtifactAbort(
            "non-DIA preview requires one exact bound artifact context",
            "multischolar_invalid_prot_nondia_artifact_context"
        )
    }
    descriptor <- findArtifactWorkflowDescriptor(
        context$getIdentity(),
        artifactWorkflowDescriptorCatalogue()
    )
    descriptor <- protNonDiaReadthroughDescriptor(descriptor$descriptor_id)
    adapter <- protNonDiaArtifactReadthroughAdapter(descriptor)
    if (is.null(import_ref)) {
        evidence <- collectProtNonDiaResumeEvidence(
            context,
            descriptor,
            resource_policy
        )
        protNonDiaValidateResumeContract(evidence)
        store <- evidence$store
        import_ref <- evidence$import$refs$canonical_data
    } else {
        if (!artifactStageIdentityMatches(context$getIdentity(), descriptor)) {
            artifactStageReadthroughAbort(
                adapter,
                "invalid_context",
                "non-DIA preview requires its exact bound artifact context"
            )
        }
        identity <- context$getIdentity()
        store <- newArtifactStore(context$getPaths(), identity$project_id)
    }
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
