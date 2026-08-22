# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

.PROT_DIA_STAGE_WORKER_SCHEMA <- "multischolar.prot_dia_stage_worker"
.PROT_DIA_STAGE_WORKER_VERSION <- 1L
.PROT_DIA_PENDING_STAGE_SCHEMA <- "multischolar.prot_dia_pending_stage"
.PROT_DIA_PENDING_STAGE_VERSION <- 1L

protDiaArtifactWorkerFailure <- function(requested_stage, current_stage) {
    if (identical(requested_stage, "timeout_before_import") &&
        identical(current_stage, "before_import")) {
        Sys.sleep(60)
    }
    if (identical(requested_stage, current_stage)) {
        protDiaArtifactAbort(
            sprintf("injected DIA-NN artifact worker failure at %s", current_stage),
            "multischolar_injected_prot_dia_worker_failure",
            failure_stage = current_stage
        )
    }
    invisible(NULL)
}

protDiaArtifactWorkerSpec <- function(
    mode,
    stage,
    source_path = NULL,
    use_precursor_norm = TRUE,
    sanitize_names = FALSE,
    failure_stage = NULL
) {
    spec <- list(
        schema = .PROT_DIA_STAGE_WORKER_SCHEMA,
        schema_version = .PROT_DIA_STAGE_WORKER_VERSION,
        mode = mode,
        stage = stage,
        source_path = source_path,
        use_precursor_norm = use_precursor_norm,
        sanitize_names = sanitize_names,
        failure_stage = failure_stage
    )
    artifactResourceDataOnly(spec, "DIA-NN stage worker specification")
    spec
}

validateProtDiaArtifactWorkerSpec <- function(spec, mode) {
    valid <- is.list(spec) &&
        identical(spec$schema, .PROT_DIA_STAGE_WORKER_SCHEMA) &&
        identical(as.integer(spec$schema_version), .PROT_DIA_STAGE_WORKER_VERSION) &&
        identical(spec$mode, mode) &&
        is.list(spec$stage) &&
        identical(spec$stage$stage_id, "import") &&
        is.list(spec$stage$store) &&
        is.list(spec$stage$identity) &&
        workflowCapabilityScalarString(spec$stage$run_id) &&
        workflowCapabilityScalarString(spec$stage$action_id) &&
        workflowCapabilityScalarString(spec$stage$generation_id)
    if (!isTRUE(valid)) {
        protDiaArtifactAbort(
            "DIA-NN artifact worker specification is malformed",
            "multischolar_invalid_prot_dia_worker_spec"
        )
    }
    artifactResourceDataOnly(spec, "DIA-NN stage worker specification")
    spec
}

sanitizeProtDiaArtifactImport <- function(data_import_result) {
    run_col <- data_import_result$column_mapping$run_col
    if (is.null(run_col) || !run_col %in% names(data_import_result$data)) {
        return(data_import_result)
    }
    original_runs <- unique(as.character(data_import_result$data[[run_col]]))
    clean_runs <- janitor::make_clean_names(original_runs)
    mapping <- stats::setNames(clean_runs, original_runs)
    data_import_result$data[[run_col]] <- unname(mapping[
        as.character(data_import_result$data[[run_col]])
    ])
    data_import_result
}

runProtDiaArtifactWriterWorker <- function(spec) {
    spec <- validateProtDiaArtifactWorkerSpec(spec, "writer")
    protDiaArtifactWorkerFailure(spec$failure_stage, "before_import")
    imported <- importDIANNData(
        spec$source_path,
        use_precursor_norm = spec$use_precursor_norm
    )
    if (isTRUE(spec$sanitize_names)) {
        imported <- sanitizeProtDiaArtifactImport(imported)
    }
    protDiaArtifactWorkerFailure(spec$failure_stage, "after_import")
    exact_digests <- list(
        canonical_data = artifactExactHydrationDigest(imported$data)
    )
    refs <- protDiaArtifactWriteStageRefs(
        spec$stage,
        list(canonical_data = imported$data),
        verification = "deferred_exact"
    )
    protDiaArtifactWorkerFailure(spec$failure_stage, "after_write")
    if (identical(spec$failure_stage, "hard_exit_after_write")) {
        quit(save = "no", status = 70L, runLast = FALSE)
    }
    source <- protDiaArtifactSourceMetadata(spec$source_path)
    result <- list(
        ok = TRUE,
        mode = "writer",
        worker_pid = as.integer(Sys.getpid()),
        refs = lapply(refs, unclass),
        exact_digests = exact_digests,
        data_type = imported$data_type,
        column_mapping = imported$column_mapping,
        source = source
    )
    protDiaArtifactWorkerFailure(spec$failure_stage, "before_result")
    artifactResourceDataOnly(result, "DIA-NN stage writer result")
    result
}

runProtDiaArtifactVerifierWorker <- function(spec) {
    spec <- validateProtDiaArtifactWorkerSpec(spec, "verifier")
    protDiaArtifactWorkerFailure(spec$failure_stage, "before_verify")
    proofs <- lapply(
        spec$stage$refs,
        \(ref) artifactStoreVerifyExactRef(spec$stage$store, ref)
    )
    names(proofs) <- names(spec$stage$refs)
    protDiaArtifactWorkerFailure(spec$failure_stage, "after_verify")
    result <- list(
        ok = TRUE,
        mode = "verifier",
        worker_pid = as.integer(Sys.getpid()),
        proofs = proofs
    )
    artifactResourceDataOnly(result, "DIA-NN stage verifier result")
    result
}

runProtDiaArtifactStageWorkerCommand <- function(spec_path, result_path, mode) {
    result <- tryCatch(
        {
            spec <- readRDS(spec_path)
            switch(mode,
                writer = runProtDiaArtifactWriterWorker(spec),
                verifier = runProtDiaArtifactVerifierWorker(spec),
                protDiaArtifactAbort(
                    "DIA-NN artifact worker mode is unsupported",
                    "multischolar_invalid_prot_dia_worker_mode"
                )
            )
        },
        error = \(error) {
            list(
                ok = FALSE,
                mode = mode,
                worker_pid = as.integer(Sys.getpid()),
                error_class = class(error),
                error_message = conditionMessage(error)
            )
        }
    )
    saveRDS(result, result_path, version = 3L)
    if (!isTRUE(result$ok)) stop(result$error_message, call. = FALSE)
    invisible(TRUE)
}

protDiaArtifactWorkerExpression <- function(source_tree) {
    load_package <- if (isTRUE(source_tree)) {
        paste0(
            "pkgload::load_all(args[[1L]], quiet = TRUE, ",
            "export_all = FALSE)"
        )
    } else {
        paste0(
            "worker_env <- new.env(parent = baseenv()); ",
            "invisible(lazyLoad(file.path(args[[1L]], 'R', 'MultiScholaR'), ",
            "envir = worker_env)); ",
            "worker_env$log_info <- logger::log_info; ",
            "worker_env$log_error <- logger::log_error; ",
            "suppressMessages(attach(worker_env, name = 'multischolar-worker'))"
        )
    }
    paste(
        "args <- commandArgs(trailingOnly = TRUE)",
        load_package,
        paste0(
            if (isTRUE(source_tree)) {
                "MultiScholaR:::runProtDiaArtifactStageWorkerCommand("
            } else {
                "worker_env$runProtDiaArtifactStageWorkerCommand("
            },
            "args[[2L]], args[[3L]], args[[4L]])"
        ),
        sep = "; "
    )
}

runProtDiaArtifactStageProcess <- function(spec, timeout_ms = 600000L) {
    if (!requireNamespace("processx", quietly = TRUE)) {
        protDiaArtifactAbort(
            "DIA-NN artifact stage workers require package 'processx'",
            "multischolar_missing_prot_dia_worker_dependency"
        )
    }
    spec <- validateProtDiaArtifactWorkerSpec(spec, spec$mode)
    spec_path <- tempfile("prot-dia-worker-spec-", fileext = ".rds")
    result_path <- tempfile("prot-dia-worker-result-", fileext = ".rds")
    on.exit(unlink(c(spec_path, result_path), force = TRUE), add = TRUE)
    saveRDS(spec, spec_path, version = 3L)
    namespace_path <- getNamespaceInfo(asNamespace("MultiScholaR"), "path")
    source_tree <- file.exists(file.path(
        namespace_path,
        "R",
        "mod_prot_artifact_stage_worker_helpers.R"
    ))
    process <- processx::run(
        command = file.path(R.home("bin"), "Rscript"),
        args = c(
            "--vanilla",
            "-e",
            protDiaArtifactWorkerExpression(source_tree),
            namespace_path,
            spec_path,
            result_path,
            spec$mode
        ),
        error_on_status = FALSE,
        timeout = as.numeric(timeout_ms) / 1000,
        echo = FALSE
    )
    if (!file.exists(result_path)) {
        protDiaArtifactAbort(
            sprintf("DIA-NN artifact %s worker exited without a result", spec$mode),
            "multischolar_prot_dia_worker_process_failed",
            worker_status = as.integer(process$status)
        )
    }
    result <- readRDS(result_path)
    artifactResourceDataOnly(result, "DIA-NN stage worker result")
    if (!isTRUE(result$ok) || !identical(as.integer(process$status), 0L)) {
        protDiaArtifactAbort(
            sprintf("DIA-NN artifact %s worker failed", spec$mode),
            "multischolar_prot_dia_worker_process_failed",
            worker_status = as.integer(process$status),
            worker_error_class = result$error_class %||% character()
        )
    }
    result
}

protDiaArtifactPreparedWorkerStage <- function(context, resource_policy = NULL) {
    stage <- protDiaArtifactStageResources(context, "import", resource_policy)
    closeProjectRegistry(stage$session)
    stage$session <- NULL
    stage$resource_policy <- normalizeProjectRegistryPolicy(resource_policy)
    artifactResourceDataOnly(stage, "DIA-NN pending import stage")
    stage
}

protDiaArtifactWorkerParameters <- function(writer, use_precursor_norm,
                                             sanitize_names) {
    list(
        column_mapping = writer$column_mapping,
        column_mapping_serialized = artifactWorkflowStateSerializeMetadata(
            writer$column_mapping,
            "DIA-NN import column mapping"
        ),
        readthrough_contract_version = 1L,
        use_precursor_norm = isTRUE(use_precursor_norm),
        sanitize_names = isTRUE(sanitize_names)
    )
}

newProtDiaArtifactPendingStage <- function(stage, writer, verifier,
                                           use_precursor_norm, sanitize_names) {
    stage$refs <- writer$refs
    stage$exact_digests <- writer$exact_digests
    pending <- list(
        schema = .PROT_DIA_PENDING_STAGE_SCHEMA,
        schema_version = .PROT_DIA_PENDING_STAGE_VERSION,
        stage = stage,
        parameters = protDiaArtifactWorkerParameters(
            writer,
            use_precursor_norm,
            sanitize_names
        ),
        source = writer$source,
        proofs = verifier$proofs,
        process_evidence = list(
            parent_pid = as.integer(Sys.getpid()),
            writer_pid = writer$worker_pid,
            verifier_pid = verifier$worker_pid
        )
    )
    validateProtDiaArtifactPendingStage(pending)
}

validateProtDiaArtifactPendingStage <- function(pending) {
    valid <- is.list(pending) &&
        identical(pending$schema, .PROT_DIA_PENDING_STAGE_SCHEMA) &&
        identical(
            as.integer(pending$schema_version),
            .PROT_DIA_PENDING_STAGE_VERSION
        ) &&
        is.list(pending$stage) &&
        identical(pending$stage$stage_id, "import") &&
        is.list(pending$stage$refs) &&
        identical(names(pending$stage$refs), names(pending$proofs)) &&
        length(pending$stage$refs) > 0L
    if (isTRUE(valid)) {
        valid <- all(vapply(names(pending$stage$refs), \(role) {
            ref <- artifactStoreNormalizeRef(pending$stage$refs[[role]])
            proof <- pending$proofs[[role]]
            identical(proof$artifact_id, ref$artifact_id) &&
                identical(
                    proof$semantic_digest,
                    ref$hash_policy$semantic$digest
                ) &&
                identical(proof$byte_digest, ref$hash_policy$byte$digest) &&
                identical(
                    proof$hydration_digest,
                    pending$stage$exact_digests[[role]]
                )
        }, logical(1)))
    }
    if (!isTRUE(valid)) {
        protDiaArtifactAbort(
            "DIA-NN pending artifact stage is malformed or unverified",
            "multischolar_invalid_prot_dia_pending_stage"
        )
    }
    artifactResourceDataOnly(pending, "DIA-NN pending import stage")
    pending
}

protDiaArtifactDiscardPath <- function(store, relative_path, recursive = FALSE) {
    if (is.null(relative_path)) return(invisible(FALSE))
    path <- artifactStoreResolveFile(store, relative_path)
    if (!file.exists(path) && !dir.exists(path)) return(invisible(FALSE))
    unlink(path, recursive = recursive, force = TRUE)
    if (file.exists(path) || dir.exists(path)) {
        protDiaArtifactAbort(
            "unregistered DIA-NN artifact path could not be discarded",
            "multischolar_prot_dia_worker_cleanup_failed"
        )
    }
    invisible(TRUE)
}

discardProtDiaArtifactWorkerStage <- function(stage) {
    store <- validateArtifactStore(stage$store)
    generation_id <- artifactRefValidateId(
        stage$generation_id,
        "generation_id",
        "gen"
    )
    for (intent_path in artifactStoreIntentPaths(store)) {
        intent <- tryCatch(artifactStoreReadIntent(store, intent_path), error = identity)
        if (inherits(intent, "error") ||
            !identical(intent$artifact_ref$logical_key$generation_id, generation_id)) {
            next
        }
        protDiaArtifactDiscardPath(
            store,
            intent$temporary_paths$payload_directory,
            recursive = TRUE
        )
        protDiaArtifactDiscardPath(store, intent$temporary_paths$sidecar)
        protDiaArtifactDiscardPath(store, intent$temporary_paths$intent)
        protDiaArtifactDiscardPath(
            store,
            intent$managed_paths$payload_directory,
            recursive = TRUE
        )
        protDiaArtifactDiscardPath(store, intent$managed_paths$sidecar)
        protDiaArtifactDiscardPath(store, intent_path)
    }
    for (role in "canonical_data") {
        protDiaArtifactDiscardPath(
            store,
            artifactNormalizeRelativePath(file.path(
                store$relative_paths$tables,
                stage$stage_id,
                role,
                generation_id
            )),
            recursive = TRUE
        )
    }
    protDiaArtifactDiscardPath(
        store,
        artifactNormalizeRelativePath(file.path(
            store$relative_paths$generations,
            generation_id
        )),
        recursive = TRUE
    )
    invisible(TRUE)
}

stageProtDiaImportArtifacts <- function(
    workflow_data,
    source_path,
    format,
    use_precursor_norm = TRUE,
    sanitize_names = FALSE,
    resource_policy = NULL,
    timeout_ms = 600000L,
    writer_failure_stage = NULL,
    verifier_failure_stage = NULL
) {
    prepared <- prepareProtDiaArtifactContext(
        workflow_data,
        format = format,
        data_type = "peptide"
    )
    if (!isTRUE(prepared$enabled)) {
        return(list(
            enabled = FALSE,
            attempted = FALSE,
            ok = TRUE,
            reason = prepared$reason
        ))
    }
    stage <- protDiaArtifactPreparedWorkerStage(
        prepared$context,
        resource_policy
    )
    caller_owns_stage <- FALSE
    on.exit({
        if (!caller_owns_stage) {
            try(discardProtDiaArtifactWorkerStage(stage), silent = TRUE)
        }
    }, add = TRUE)
    writer <- runProtDiaArtifactStageProcess(protDiaArtifactWorkerSpec(
        mode = "writer",
        stage = stage,
        source_path = source_path,
        use_precursor_norm = use_precursor_norm,
        sanitize_names = sanitize_names,
        failure_stage = writer_failure_stage
    ), timeout_ms)
    verifier_stage <- stage
    verifier_stage$refs <- writer$refs
    verifier <- runProtDiaArtifactStageProcess(protDiaArtifactWorkerSpec(
        mode = "verifier",
        stage = verifier_stage,
        failure_stage = verifier_failure_stage
    ), timeout_ms)
    pending <- newProtDiaArtifactPendingStage(
        stage,
        writer,
        verifier,
        use_precursor_norm,
        sanitize_names
    )
    tables <- readProtDiaArtifactStage(stage$store, pending$stage)
    imported <- list(
        data = tables$canonical_data,
        data_type = writer$data_type,
        column_mapping = writer$column_mapping
    )
    caller_owns_stage <- TRUE
    list(
        enabled = TRUE,
        attempted = TRUE,
        ok = TRUE,
        result = imported,
        pending_stage = pending,
        process_evidence = pending$process_evidence
    )
}

stageProtDiaImportArtifactsSafely <- function(
    workflow_data,
    source_path,
    format,
    use_precursor_norm = TRUE,
    sanitize_names = FALSE,
    resource_policy = NULL,
    timeout_ms = 600000L,
    writer_failure_stage = NULL,
    verifier_failure_stage = NULL,
    log_warn = logger::log_warn
) {
    result <- tryCatch(
        stageProtDiaImportArtifacts(
            workflow_data = workflow_data,
            source_path = source_path,
            format = format,
            use_precursor_norm = use_precursor_norm,
            sanitize_names = sanitize_names,
            resource_policy = resource_policy,
            timeout_ms = timeout_ms,
            writer_failure_stage = writer_failure_stage,
            verifier_failure_stage = verifier_failure_stage
        ),
        error = \(error) {
            log_warn(paste(
                "DIA-NN artifact import worker failed; using the memory importer:",
                conditionMessage(error)
            ))
            list(
                enabled = TRUE,
                attempted = TRUE,
                ok = FALSE,
                error_class = class(error),
                error_message = conditionMessage(error)
            )
        }
    )
    recordProtDiaArtifactResult(workflow_data, "import_worker", result)
    result
}

publishProtDiaArtifactWorkerStage <- function(workflow_data, pending_stage) {
    pending <- validateProtDiaArtifactPendingStage(pending_stage)
    context <- workflow_data$workflow_context
    if (!inherits(context, "WorkflowContext") || !context$isBound() ||
        !identical(context$getIdentity(), pending$stage$identity)) {
        protDiaArtifactAbort(
            "DIA-NN pending stage belongs to another workflow context",
            "multischolar_invalid_prot_dia_pending_stage_context"
        )
    }
    registry <- projectRegistryForContext(
        context,
        pending$stage$resource_policy
    )
    session <- initializeProjectRegistry(registry)
    on.exit(closeProjectRegistry(session), add = TRUE)
    artifactWorkflowStateEnsureWorkflow(session, pending$stage$identity)
    stage <- pending$stage
    stage$session <- session
    protDiaArtifactRegisterStage(
        stage,
        stage$refs,
        pending$parameters,
        pending$source,
        deferred_commit = FALSE
    )
    list(
        enabled = TRUE,
        ok = TRUE,
        stage_id = stage$stage_id,
        run_id = stage$run_id,
        action_id = stage$action_id,
        generation_id = stage$generation_id,
        refs = stage$refs,
        committed = TRUE,
        process_evidence = pending$process_evidence
    )
}
