# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

#' Collect coordinator-owned non-DIA design tables
#'
#' @param workflow_data Mutable proteomics workflow state.
#' @param state_object Completed memory-backed protein S4 checkpoint.
#'
#' @return A named list of rectangular design-stage tables.
#' @noRd
protNonDiaDesignArtifactTables <- function(workflow_data, state_object) {
    args <- if (isS4(state_object) && "args" %in% methods::slotNames(state_object)) {
        methods::slot(state_object, "args")
    } else {
        workflow_data$config_list
    }
    contrasts <- workflow_data$contrasts_tbl
    if (is.character(contrasts)) {
        contrasts <- data.frame(contrasts = contrasts, stringsAsFactors = FALSE)
    }
    list(
        cleaned_data = workflow_data$data_cln,
        design_matrix = workflow_data$design_matrix,
        contrasts = artifactStageOptionalTable(
            contrasts,
            "contrasts",
            protNonDiaArtifactAbort
        ),
        args = artifactStageMetadataTable(args),
        annotations = artifactStageOptionalTable(
            workflow_data$uniprot_dat_cln,
            "annotations",
            protNonDiaArtifactAbort
        ),
        sequences = artifactStageOptionalTable(
            workflow_data$aa_seq_tbl_final,
            "sequences",
            protNonDiaArtifactAbort
        )
    )
}

#' Resolve committed parent import lineage
#'
#' @param workflow_data Mutable proteomics workflow state.
#'
#' @return Import run, generation, artifact, and semantic digest identifiers.
#' @noRd
protNonDiaArtifactImportParent <- function(workflow_data) {
    result <- workflow_data$artifact_stage_results$import
    ref <- if (is.list(result)) result$refs$canonical_data else NULL
    valid <- is.list(result) && isTRUE(result$enabled) && isTRUE(result$ok) &&
        isTRUE(result$committed) && workflowCapabilityScalarString(result$run_id) &&
        workflowCapabilityScalarString(result$generation_id) && is.list(ref) &&
        workflowCapabilityScalarString(ref$artifact_id) &&
        workflowCapabilityScalarString(ref$hash_policy$semantic$digest)
    if (!isTRUE(valid)) {
        return(list(
            run_id = NULL,
            generation_id = NULL,
            artifact_id = NULL,
            semantic_digest = NULL
        ))
    }
    list(
        run_id = result$run_id,
        generation_id = result$generation_id,
        artifact_id = ref$artifact_id,
        semantic_digest = ref$hash_policy$semantic$digest
    )
}

#' Validate the completed non-DIA protein memory checkpoint
#'
#' @param workflow_data Mutable proteomics workflow state.
#' @param state_name Requested workflow state name.
#' @param prepared Prepared exact-tuple artifact context.
#'
#' @return The memory state manager and validated S4 checkpoint.
#' @noRd
protNonDiaArtifactMemoryCheckpoint <- function(workflow_data, state_name, prepared) {
    if (!identical(state_name, "protein_s4_initial")) {
        protNonDiaArtifactAbort(
            "non-DIA design dual-write requires 'protein_s4_initial'",
            "multischolar_unsupported_proteomics_codec_role",
            state_role = state_name
        )
    }
    manager <- workflow_data$state_manager
    state_object <- manager$getState(state_name)
    if (is.null(state_object) || !isS4(state_object)) {
        protNonDiaArtifactAbort(
            "non-DIA design dual-write requires the completed memory checkpoint",
            "multischolar_missing_prot_nondia_memory_checkpoint"
        )
    }
    artifactValidateProteomicsNonDiaProteinState(
        state_object,
        artifactProteomicsNonDiaCodecRole(
            prepared$descriptor$descriptor_id,
            state_name,
            prepared$descriptor$descriptor_version
        ),
        state_name
    )
    list(manager = manager, state_object = state_object)
}

#' Write a deferred non-DIA design artifact generation
#'
#' @param prepared Prepared exact-tuple artifact context.
#' @param workflow_data Mutable proteomics workflow state.
#' @param state_name Completed memory state name.
#' @param checkpoint Validated memory checkpoint.
#' @param failure_injector Optional artifact failure injector used by tests.
#'
#' @return A validated but uncommitted artifact stage record.
#' @noRd
protNonDiaArtifactWriteDesignStage <- function(
    prepared,
    workflow_data,
    state_name,
    checkpoint,
    failure_injector
) {
    parent_import <- protNonDiaArtifactImportParent(workflow_data)
    stage <- writeArtifactStage(
        prepared$context,
        prepared$descriptor,
        stage_id = "design",
        tables = protNonDiaDesignArtifactTables(
            workflow_data,
            checkpoint$state_object
        ),
        parameters = list(
            capability_id = prepared$descriptor$descriptor_id,
            state_name = state_name,
            workflow_type = workflowStateType(checkpoint$manager),
            readthrough_contract_version = 1L,
            parent_import_run_id = parent_import$run_id,
            parent_import_generation_id = parent_import$generation_id,
            parent_import_artifact_id = parent_import$artifact_id,
            parent_import_semantic_digest = parent_import$semantic_digest,
            contrasts_kind = artifactStageContrastsKind(
                workflow_data$contrasts_tbl
            ),
            annotation_available = !is.null(workflow_data$uniprot_dat_cln),
            sequence_data_available = !is.null(workflow_data$aa_seq_tbl_final)
        ),
        deferred_commit = TRUE,
        failure_injector = failure_injector,
        abort_fn = protNonDiaArtifactAbort
    )
    stage$readthrough_contract_version <- 1L
    stage$parent_import <- parent_import
    stage
}

#' Build audit metadata for a non-DIA design checkpoint
#'
#' @param stage Validated design artifact stage.
#'
#' @return Audit metadata linking the S4 checkpoint and stage artifacts.
#' @noRd
protNonDiaArtifactStateAudit <- function(stage) {
    list(
        record_id = artifactOpaqueId("audit"),
        provenance_status = "artifact_dual_write_validated",
        stage_id = stage$stage_id,
        run_id = stage$run_id,
        action_id = stage$action_id,
        readthrough_contract_version = stage$readthrough_contract_version,
        parent_import_run_id = stage$parent_import$run_id,
        parent_import_generation_id = stage$parent_import$generation_id,
        parent_import_artifact_id = stage$parent_import$artifact_id,
        parent_import_semantic_digest = stage$parent_import$semantic_digest,
        stage_artifact_refs = stage$refs
    )
}

#' Construct a descriptor-bound artifact state manager
#'
#' @param prepared Prepared exact-tuple artifact context.
#' @param manager_factory Workflow state manager constructor.
#'
#' @return A workflow state manager bound to the non-DIA descriptor.
#' @noRd
protNonDiaArtifactNewDesignStateManager <- function(prepared, manager_factory) {
    manager_factory(
        workflow_context = prepared$context,
        workflow_descriptor = prepared$descriptor,
        descriptor_catalogue = artifactWorkflowDescriptorCatalogue(),
        codec_catalogue = artifactS4CodecCatalogue()
    )
}

#' Persist and independently hydrate a non-DIA protein S4 checkpoint
#'
#' @param manager Descriptor-bound artifact state manager.
#' @param workflow_data Mutable proteomics workflow state.
#' @param state_name Completed memory state name.
#' @param state_object Completed memory-backed protein S4 checkpoint.
#' @param stage Validated design artifact stage.
#' @param prepared Prepared exact-tuple artifact context.
#' @param save_state_fn Optional state writer injected by tests.
#'
#' @return State manifest and metadata from independent hydration.
#' @noRd
protNonDiaArtifactSaveDesignState <- function(
    manager,
    workflow_data,
    state_name,
    state_object,
    stage,
    prepared,
    save_state_fn
) {
    spec <- protNonDiaArtifactImportSpec(workflow_data$data_format)
    manager$setWorkflowType(spec$workflow_type)
    audit <- protNonDiaArtifactStateAudit(stage)
    if (is.null(save_state_fn)) {
        manager$saveState(
            state_name = state_name,
            s4_data_object = state_object,
            config_object = workflow_data$config_list,
            description = sprintf(
                "%s non-DIA protein design dual-write checkpoint.",
                spec$workflow_type
            ),
            audit_metadata = audit
        )
    } else {
        save_state_fn(manager, state_name, state_object, audit)
    }
    hydrated <- manager$getState(state_name)
    artifactValidateProteomicsNonDiaProteinState(
        hydrated,
        artifactProteomicsNonDiaCodecRole(
            prepared$descriptor$descriptor_id,
            state_name,
            prepared$descriptor$descriptor_version
        ),
        state_name
    )
    if (!identical(hydrated, state_object)) {
        protNonDiaArtifactAbort(
            "non-DIA design artifact differs from its memory checkpoint",
            "multischolar_inexact_prot_nondia_artifact_hydration"
        )
    }
    list(
        state_manifest = manager$exportState(),
        state_metadata = manager$getStateMetadata(state_name)
    )
}

#' Fail a deferred non-DIA design stage
#'
#' @param context Bound workflow context.
#' @param stage Deferred design artifact stage.
#' @param state_error Original state persistence or hydration error.
#'
#' @return This function does not return; it re-signals `state_error`.
#' @noRd
protNonDiaArtifactFailDesignStage <- function(context, stage, state_error) {
    try(
        artifactStageUpdateStatus(
            context,
            stage,
            completed = FALSE,
            abort_fn = protNonDiaArtifactAbort
        ),
        silent = TRUE
    )
    stop(state_error)
}

#' Digest one non-DIA scientific checkpoint exactly
#'
#' @param object Validated protein S4 checkpoint.
#' @param descriptor Exact non-DIA workflow descriptor.
#' @param state_name Exact state role.
#'
#' @return A SHA-256 digest string.
#' @noRd
protNonDiaDesignStateDigest <- function(object, descriptor, state_name) {
    slot_names <- methods::slotNames(object)
    slot_digests <- vapply(slot_names, \(slot_name) {
        value <- methods::slot(object, slot_name)
        if (is.data.frame(value)) return(artifactExactHydrationDigest(value))
        digest::digest(
            value,
            algo = "sha256",
            serialize = TRUE,
            ascii = FALSE,
            serializeVersion = 3L
        )
    }, character(1))
    artifactSemanticDigest(list(
        descriptor_id = descriptor$descriptor_id,
        descriptor_version = descriptor$descriptor_version,
        state_name = state_name,
        class_name = class(object)[[1L]],
        slot_names = slot_names,
        slot_classes = vapply(slot_names, \(slot_name) {
            class(methods::slot(object, slot_name))[[1L]]
        }, character(1)),
        slot_digests = slot_digests
    ))
}

#' Write one non-DIA design checkpoint without hydration
#'
#' @param workflow_data Mutable proteomics workflow state.
#' @param state_name Completed memory state name.
#'
#' @return A payload-free pending design result.
#' @noRd
runProtNonDiaDesignWriter <- function(workflow_data, state_name) {
    prepared <- prepareProtNonDiaArtifactContext(workflow_data)
    checkpoint <- protNonDiaArtifactMemoryCheckpoint(
        workflow_data,
        state_name,
        prepared
    )
    stage <- protNonDiaArtifactWriteDesignStage(
        prepared,
        workflow_data,
        state_name,
        checkpoint,
        failure_injector = NULL
    )
    manager <- protNonDiaArtifactNewDesignStateManager(
        prepared,
        newWorkflowState
    )
    on.exit(manager$close(), add = TRUE)
    spec <- protNonDiaArtifactImportSpec(workflow_data$data_format)
    manager$setWorkflowType(spec$workflow_type)
    audit <- protNonDiaArtifactStateAudit(stage)
    state_error <- tryCatch(
        {
            manager$saveState(
                state_name = state_name,
                s4_data_object = checkpoint$state_object,
                config_object = workflow_data$config_list,
                description = sprintf(
                    "%s non-DIA protein design worker checkpoint.",
                    spec$workflow_type
                ),
                audit_metadata = audit
            )
            NULL
        },
        error = identity
    )
    if (inherits(state_error, "error")) {
        protNonDiaArtifactFailDesignStage(
            prepared$context,
            stage,
            state_error
        )
    }
    result <- list(
        ok = TRUE,
        worker_pid = as.integer(Sys.getpid()),
        stage = stage,
        state_manifest = manager$exportState(),
        state_metadata = manager$getStateMetadata(state_name),
        expected_state_digest = protNonDiaDesignStateDigest(
            checkpoint$state_object,
            prepared$descriptor,
            state_name
        ),
        complete_payload_returned = FALSE
    )
    artifactResourceDataOnly(result, "non-DIA design writer result")
    result
}

#' Independently hydrate one non-DIA design checkpoint
#'
#' @param workflow_data Mutable proteomics workflow state.
#' @param state_name Expected current state name.
#' @param expected_digest Exact writer-side scientific digest.
#'
#' @return A payload-free hydration proof.
#' @noRd
runProtNonDiaDesignVerifier <- function(
    workflow_data,
    state_name,
    expected_digest
) {
    prepared <- prepareProtNonDiaArtifactContext(workflow_data)
    manager <- protNonDiaArtifactNewDesignStateManager(
        prepared,
        newWorkflowState
    )
    on.exit(manager$close(), add = TRUE)
    hydrated <- manager$getState(state_name)
    artifactValidateProteomicsNonDiaProteinState(
        hydrated,
        artifactProteomicsNonDiaCodecRole(
            prepared$descriptor$descriptor_id,
            state_name,
            prepared$descriptor$descriptor_version
        ),
        state_name
    )
    hydrated_digest <- protNonDiaDesignStateDigest(
        hydrated,
        prepared$descriptor,
        state_name
    )
    if (!identical(hydrated_digest, expected_digest)) {
        protNonDiaArtifactAbort(
            sprintf(
                paste(
                    "non-DIA design verifier differs from the memory",
                    "checkpoint (%s != %s)"
                ),
                hydrated_digest,
                expected_digest
            ),
            "multischolar_inexact_prot_nondia_artifact_hydration"
        )
    }
    manager$releaseCache()
    result <- list(
        ok = TRUE,
        worker_pid = as.integer(Sys.getpid()),
        state_name = state_name,
        hydrated_digest = hydrated_digest,
        valid_s4 = identical(methods::validObject(hydrated, test = TRUE), TRUE),
        complete_payload_returned = FALSE
    )
    artifactResourceDataOnly(result, "non-DIA design verifier result")
    result
}

#' Run one non-DIA design operation in a fork child
#'
#' @param operation Zero-argument child operation.
#' @param owner Human-readable process owner.
#'
#' @return A payload-free child result.
#' @noRd
runProtNonDiaDesignFork <- function(
    operation,
    owner,
    timeout_ms = 900000L
) {
    job <- parallel::mcparallel(
        tryCatch(
            operation(),
            error = \(error) list(
                ok = FALSE,
                worker_pid = as.integer(Sys.getpid()),
                error_class = class(error),
                error_message = conditionMessage(error)
            )
        ),
        mc.set.seed = FALSE,
        silent = TRUE
    )
    deadline <- proc.time()[["elapsed"]] + as.numeric(timeout_ms) / 1000
    result <- NULL
    repeat {
        collected <- parallel::mccollect(job, wait = FALSE)
        if (!is.null(collected)) {
            result <- unname(collected)[[1L]]
            break
        }
        if (proc.time()[["elapsed"]] >= deadline) {
            tools::pskill(job$pid, signal = 15L)
            parallel::mccollect(job, wait = TRUE)
            protNonDiaArtifactAbort(
                paste(owner, "exceeded its timeout"),
                "multischolar_prot_nondia_design_worker_failed"
            )
        }
        Sys.sleep(0.05)
    }
    if (inherits(result, "try-error") || !is.list(result) ||
        !isTRUE(result$ok)) {
        error_class <- if (is.list(result)) {
            result$error_class %||% character()
        } else {
            class(result)
        }
        error_message <- if (is.list(result)) {
            result$error_message %||% as.character(result)
        } else {
            as.character(result)
        }
        protNonDiaArtifactAbort(
            paste(owner, "failed"),
            unique(c(
                error_class,
                "multischolar_prot_nondia_design_worker_failed"
            )),
            worker_error_class = error_class,
            worker_error_message = error_message
        )
    }
    artifactResourceDataOnly(result, owner)
    result
}

#' Roll back an unverified non-DIA design state
#'
#' @param workflow_data Mutable proteomics workflow state.
#' @param stage Pending design stage.
#'
#' @return `FALSE`, invisibly.
#' @noRd
rollbackProtNonDiaDesignVerification <- function(workflow_data, stage) {
    try(
        artifactStageUpdateStatus(
            workflow_data$workflow_context,
            stage,
            completed = FALSE,
            abort_fn = protNonDiaArtifactAbort
        ),
        silent = TRUE
    )
    descriptor <- protNonDiaReadthroughDescriptor(
        stage$descriptor$descriptor_id
    )
    manager <- tryCatch(
        newProtNonDiaResumeStateManager(
            workflow_data$workflow_context,
            descriptor
        ),
        error = identity
    )
    if (inherits(manager, "ArtifactWorkflowState")) {
        on.exit(manager$close(), add = TRUE)
        if (isTRUE(manager$hasState("initial"))) {
            try(manager$revertToState("initial"), silent = TRUE)
        }
    }
    invisible(FALSE)
}

#' Persist and verify one non-DIA design in separate processes
#'
#' @param workflow_data Mutable proteomics workflow state.
#' @param state_name Completed memory state name.
#'
#' @return A committed payload-free design result.
#' @noRd
persistProtNonDiaDesignArtifactsProcessCore <- function(
    workflow_data,
    state_name
) {
    writer <- runProtNonDiaDesignFork(
        \() runProtNonDiaDesignWriter(workflow_data, state_name),
        "non-DIA design writer"
    )
    verifier <- tryCatch(
        runProtNonDiaDesignFork(
            \() runProtNonDiaDesignVerifier(
                workflow_data,
                state_name,
                writer$expected_state_digest
            ),
            "non-DIA design verifier"
        ),
        error = identity
    )
    if (inherits(verifier, "error")) {
        rollbackProtNonDiaDesignVerification(workflow_data, writer$stage)
        stop(verifier)
    }
    artifactStageUpdateStatus(
        workflow_data$workflow_context,
        writer$stage,
        completed = TRUE,
        abort_fn = protNonDiaArtifactAbort
    )
    writer$stage$committed <- TRUE
    c(
        list(enabled = TRUE, ok = TRUE),
        writer$stage,
        list(
            state_manifest = writer$state_manifest,
            state_metadata = writer$state_metadata,
            process_evidence = list(
                writer_pid = writer$worker_pid,
                verifier_pid = verifier$worker_pid,
                distinct_workers = !identical(
                    writer$worker_pid,
                    verifier$worker_pid
                ),
                valid_s4 = verifier$valid_s4,
                exact_digest = verifier$hydrated_digest,
                complete_payload_returned = FALSE
            )
        )
    )
}

#' Execute non-DIA design dual-write transaction mechanics
#'
#' @param workflow_data Mutable proteomics workflow state.
#' @param state_name Completed memory state name.
#' @param failure_injector Optional artifact failure injector used by tests.
#' @param manager_factory Workflow state manager constructor.
#' @param save_state_fn Optional state writer injected by tests.
#'
#' @return A committed design stage and independently hydrated state metadata.
#' @noRd
persistProtNonDiaDesignArtifactsCore <- function(
    workflow_data,
    state_name,
    failure_injector = NULL,
    manager_factory = newWorkflowState,
    save_state_fn = NULL
) {
    prepared <- prepareProtNonDiaArtifactContext(workflow_data)
    if (!isTRUE(prepared$enabled)) {
        return(list(
            enabled = FALSE,
            ok = TRUE,
            stage_id = "design",
            reason = prepared$reason
        ))
    }
    checkpoint <- protNonDiaArtifactMemoryCheckpoint(
        workflow_data,
        state_name,
        prepared
    )
    stage <- protNonDiaArtifactWriteDesignStage(
        prepared,
        workflow_data,
        state_name,
        checkpoint,
        failure_injector
    )
    manager <- NULL
    state_result <- tryCatch(
        {
            manager <- protNonDiaArtifactNewDesignStateManager(
                prepared,
                manager_factory
            )
            on.exit({
                if (!is.null(manager)) manager$close()
            }, add = TRUE)
            protNonDiaArtifactSaveDesignState(
                manager,
                workflow_data,
                state_name,
                checkpoint$state_object,
                stage,
                prepared,
                save_state_fn
            )
        },
        error = \(error) error
    )
    if (inherits(state_result, "error")) {
        if (!is.null(manager)) {
            try(manager$close(), silent = TRUE)
            manager <- NULL
        }
        protNonDiaArtifactFailDesignStage(
            prepared$context,
            stage,
            state_result
        )
    }
    manager$close()
    manager <- NULL
    artifactStageUpdateStatus(
        prepared$context,
        stage,
        completed = TRUE,
        failure_injector = failure_injector,
        abort_fn = protNonDiaArtifactAbort
    )
    stage$committed <- TRUE
    c(list(enabled = TRUE, ok = TRUE), stage, state_result)
}

#' Safely dual-write a completed non-DIA proteomics design
#'
#' @param workflow_data Mutable proteomics workflow state.
#' @param state_name Completed memory state name.
#' @param failure_injector Optional artifact failure injector used by tests.
#' @param manager_factory Workflow state manager constructor.
#' @param save_state_fn Optional state writer injected by tests.
#' @param log_warn Warning logger used for additive artifact failures.
#'
#' @return The recorded design persistence result, invisibly.
#' @noRd
persistProtNonDiaDesignArtifacts <- function(
    workflow_data,
    state_name = workflowStateCurrentName(workflow_data$state_manager),
    failure_injector = NULL,
    manager_factory = newWorkflowState,
    save_state_fn = NULL,
    log_warn = logger::log_warn
) {
    process_owned <- workflow_data$data_format %in% c("maxquant", "fragpipe") &&
        identical(.Platform$OS.type, "unix") &&
        is.null(failure_injector) &&
        identical(manager_factory, newWorkflowState) &&
        is.null(save_state_fn)
    runArtifactStageSafely(
        workflow_data,
        "design",
        if (isTRUE(process_owned)) {
            \() persistProtNonDiaDesignArtifactsProcessCore(
                workflow_data,
                state_name
            )
        } else {
            \() persistProtNonDiaDesignArtifactsCore(
                workflow_data,
                state_name,
                failure_injector,
                manager_factory,
                save_state_fn
            )
        },
        "non-DIA proteomics",
        log_warn
    )
}

#' Dispatch proteomics design dual-write by exact workflow tuple
#'
#' @param workflow_data Mutable proteomics workflow state.
#' @param state_name Completed memory state name.
#' @param failure_injector Optional artifact failure injector used by tests.
#' @param manager_factory Workflow state manager constructor.
#' @param save_state_fn Optional state writer injected by tests.
#' @param log_warn Warning logger used for additive artifact failures.
#'
#' @return The recorded DIA or non-DIA design persistence result, invisibly.
#' @noRd
persistProtDesignArtifacts <- function(
    workflow_data,
    state_name = workflowStateCurrentName(workflow_data$state_manager),
    failure_injector = NULL,
    manager_factory = newWorkflowState,
    save_state_fn = NULL,
    log_warn = logger::log_warn
) {
    if (identical(workflow_data$data_format, "diann") &&
        identical(workflow_data$data_type, "peptide")) {
        return(persistProtDiaDesignArtifacts(
            workflow_data,
            state_name,
            failure_injector,
            manager_factory,
            save_state_fn,
            log_warn
        ))
    }
    persistProtNonDiaDesignArtifacts(
        workflow_data,
        state_name,
        failure_injector,
        manager_factory,
        save_state_fn,
        log_warn
    )
}
