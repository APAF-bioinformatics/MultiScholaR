# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

.PROT_DIA_DESIGN_WORKER_SCHEMA <- "multischolar.prot_dia_design_worker"
.PROT_DIA_DESIGN_WORKER_VERSION <- 1L

protDiaDeferredDesignSpec <- function(
    workflow_data,
    results,
    experiment_paths,
    experiment_label
) {
    mapping <- protDiaDeferredRunMapping(results$data_cln)
    import <- workflow_data$artifact_stage_results$import
    if (is.null(mapping) || !is.list(import) || !isTRUE(import$committed)) {
        protDiaArtifactAbort(
            "DIA-NN deferred design lacks run mapping or import evidence",
            "multischolar_invalid_prot_dia_design_worker_spec"
        )
    }
    spec <- list(
        schema = .PROT_DIA_DESIGN_WORKER_SCHEMA,
        schema_version = .PROT_DIA_DESIGN_WORKER_VERSION,
        experiment_paths = experiment_paths,
        experiment_label = experiment_label,
        storage_policy = workflow_data$workflow_context$getStoragePolicy(),
        import_result = import,
        run_mapping = mapping,
        design_matrix = results$design_matrix,
        contrasts_tbl = results$contrasts_tbl,
        config_list = results$config_list,
        uniprot_dat_cln = workflow_data$uniprot_dat_cln,
        aa_seq_tbl_final = workflow_data$aa_seq_tbl_final
    )
    artifactResourceDataOnly(spec, "DIA-NN design worker specification")
    spec
}

validateProtDiaDeferredDesignSpec <- function(spec) {
    required <- c(
        "schema", "schema_version", "experiment_paths", "experiment_label",
        "storage_policy", "import_result", "run_mapping", "design_matrix",
        "contrasts_tbl", "config_list", "uniprot_dat_cln", "aa_seq_tbl_final"
    )
    valid <- is.list(spec) && identical(names(spec), required) &&
        identical(spec$schema, .PROT_DIA_DESIGN_WORKER_SCHEMA) &&
        identical(as.integer(spec$schema_version), .PROT_DIA_DESIGN_WORKER_VERSION) &&
        is.list(spec$experiment_paths) &&
        workflowCapabilityScalarString(spec$experiment_label) &&
        is.list(spec$storage_policy) && is.list(spec$import_result) &&
        isTRUE(spec$import_result$committed) && is.data.frame(spec$run_mapping) &&
        is.data.frame(spec$design_matrix) && is.list(spec$config_list)
    if (!isTRUE(valid)) {
        protDiaArtifactAbort(
            "DIA-NN design worker specification is malformed",
            "multischolar_invalid_prot_dia_design_worker_spec"
        )
    }
    artifactResourceDataOnly(spec, "DIA-NN design worker specification")
    spec
}

protDiaDeferredDesignWorkflow <- function(spec) {
    context <- createWorkflowContext(
        spec$experiment_paths,
        "proteomics",
        spec$experiment_label,
        storage_policy = spec$storage_policy
    )
    workflow <- new.env(parent = emptyenv())
    workflow$workflow_context <- context
    workflow$state_manager <- WorkflowState$new()
    workflow$state_manager$setWorkflowType("DIA")
    workflow$artifact_stage_results <- list(import = spec$import_result)
    workflow$data_format <- "diann"
    workflow$data_type <- "peptide"
    workflow$column_mapping <- spec$import_result$column_mapping
    prepareProtDiaArtifactContext(workflow)
    workflow
}

protDiaDeferredDesignData <- function(workflow, spec) {
    ref <- spec$import_result$refs$canonical_data
    store <- newArtifactStore(
        workflow$workflow_context$getPaths(),
        workflow$workflow_context$getIdentity()$project_id
    )
    data <- protDiaArtifactReadTable(store, ref)
    mapping <- spec$run_mapping
    index <- match(as.character(data$Run), mapping$original_run)
    keep <- !is.na(index)
    data <- data[keep, , drop = FALSE]
    data$Run <- mapping$Run[index[keep]]
    tibble::as_tibble(data)
}

protDiaDeferredDesignManagerFactory <- function(..., verify_hydration_fn) {
    newWorkflowState(
        ...,
        verify_hydration_fn = artifactWorkflowStateVerifyDeferredDigest
    )
}

runProtDiaDeferredDesignWorker <- function(spec) {
    spec <- validateProtDiaDeferredDesignSpec(spec)
    workflow <- protDiaDeferredDesignWorkflow(spec)
    data <- protDiaDeferredDesignData(workflow, spec)
    workflow$data_tbl <- NULL
    workflow$data_cln <- data
    workflow$design_matrix <- spec$design_matrix
    workflow$contrasts_tbl <- spec$contrasts_tbl
    workflow$config_list <- spec$config_list
    workflow$uniprot_dat_cln <- spec$uniprot_dat_cln
    workflow$aa_seq_tbl_final <- spec$aa_seq_tbl_final
    object <- PeptideQuantitativeDataDiann(
        peptide_data = data,
        design_matrix = spec$design_matrix,
        technical_replicate_id = "tech_rep_group",
        args = spec$config_list
    )
    workflow$state_manager$saveState(
        "raw_data_s4",
        object,
        spec$config_list,
        "DIA-NN deferred design worker checkpoint."
    )
    compatibility <- list(
        data_cln = data,
        design_matrix = spec$design_matrix,
        contrasts_tbl = spec$contrasts_tbl,
        config_list = spec$config_list
    )
    persistProtDesignBuilderArtifacts(
        compatibility,
        workflow,
        spec$experiment_paths$source_dir
    )
    design <- persistProtDiaDesignArtifacts(
        workflow,
        manager_factory = protDiaDeferredDesignManagerFactory
    )
    if (!isTRUE(design$ok) || !isTRUE(design$committed)) {
        protDiaArtifactAbort(
            "DIA-NN deferred design artifact publication failed",
            "multischolar_prot_dia_design_worker_failed"
        )
    }
    result <- list(
        ok = TRUE,
        worker_pid = as.integer(Sys.getpid()),
        design_result = design,
        readthrough_proof = NULL,
        compatibility_checkpoint = NULL,
        import_refs = spec$import_result$refs,
        design_refs = design$refs,
        state_generation_id = design$state_manifest$current_generation_id,
        exact_state_digest = design$hydration_verification$expected_digest,
        complete_payload_returned = FALSE
    )
    artifactResourceDataOnly(result, "DIA-NN design worker result")
    result
}

runProtDiaDeferredDesignWorkerCommand <- function(spec_path, result_path) {
    result <- tryCatch(
        runProtDiaDeferredDesignWorker(readRDS(spec_path)),
        error = function(error) list(
            ok = FALSE,
            worker_pid = as.integer(Sys.getpid()),
            error_class = class(error),
            error_message = conditionMessage(error),
            complete_payload_returned = FALSE
        )
    )
    saveRDS(result, result_path, version = 3L)
    if (!isTRUE(result$ok)) stop(result$error_message, call. = FALSE)
    invisible(TRUE)
}

protDiaDeferredDesignWorkerExpression <- function(source_tree) {
    load_package <- if (isTRUE(source_tree)) {
        paste0(
            "pkgload::load_all(args[[1L]], quiet = TRUE, ",
            "export_all = FALSE)"
        )
    } else {
        "worker_env <- loadNamespace('MultiScholaR', lib.loc = dirname(args[[1L]]))"
    }
    call <- if (isTRUE(source_tree)) {
        "MultiScholaR:::runProtDiaDeferredDesignWorkerCommand"
    } else {
        "worker_env$runProtDiaDeferredDesignWorkerCommand"
    }
    paste(
        "args <- commandArgs(trailingOnly = TRUE)",
        load_package,
        paste0(call, "(args[[2L]], args[[3L]])"),
        sep = "; "
    )
}

protDiaDeferredDesignWorkerMode <- function() {
    requested <- getOption(
        "multischolar.prot_dia.design_worker_mode",
        if (identical(Sys.info()[["sysname"]], "Linux")) "fork" else "process"
    )
    match.arg(requested, c("fork", "process"))
}

protDiaDeferredDesignParityMode <- function() {
    requested <- getOption(
        "multischolar.prot_dia.design_parity_worker_mode",
        if (identical(Sys.info()[["sysname"]], "Linux")) "fork" else "process"
    )
    match.arg(requested, c("fork", "process"))
}

runProtDiaDeferredDesignFork <- function(spec, timeout_ms) {
    job <- parallel::mcparallel(
        tryCatch(
            runProtDiaDeferredDesignWorker(spec),
            error = function(error) list(
                ok = FALSE,
                worker_pid = as.integer(Sys.getpid()),
                error_class = class(error),
                error_message = conditionMessage(error),
                complete_payload_returned = FALSE
            )
        ),
        mc.set.seed = FALSE,
        silent = TRUE
    )
    deadline <- proc.time()[["elapsed"]] + as.numeric(timeout_ms) / 1000
    repeat {
        collected <- parallel::mccollect(job, wait = FALSE)
        if (!is.null(collected)) return(unname(collected)[[1L]])
        if (proc.time()[["elapsed"]] >= deadline) {
            tools::pskill(job$pid, signal = 15L)
            parallel::mccollect(job, wait = TRUE)
            protDiaArtifactAbort(
                "DIA-NN design worker exceeded its timeout",
                "multischolar_prot_dia_design_worker_failed"
            )
        }
        Sys.sleep(0.05)
    }
}

runProtDiaDeferredDesignProcessx <- function(spec, timeout_ms) {
    spec_path <- tempfile("prot-dia-design-spec-", fileext = ".rds")
    result_path <- tempfile("prot-dia-design-result-", fileext = ".rds")
    on.exit(unlink(c(spec_path, result_path), force = TRUE), add = TRUE)
    saveRDS(spec, spec_path, version = 3L)
    namespace_path <- getNamespaceInfo(asNamespace("MultiScholaR"), "path")
    source_tree <- file.exists(file.path(
        namespace_path,
        "R",
        "mod_prot_dia_design_worker_helpers.R"
    ))
    process <- processx::run(
        command = file.path(R.home("bin"), "Rscript"),
        args = c(
            "--vanilla", "-e", protDiaDeferredDesignWorkerExpression(source_tree),
            namespace_path, spec_path, result_path
        ),
        error_on_status = FALSE,
        timeout = as.numeric(timeout_ms) / 1000,
        echo = FALSE
    )
    if (!file.exists(result_path)) {
        protDiaArtifactAbort(
            "DIA-NN design worker exited without a result",
            "multischolar_prot_dia_design_worker_failed"
        )
    }
    list(result = readRDS(result_path), status = as.integer(process$status))
}

protDiaDeferredDesignParitySpec <- function(spec, result) {
    workflow <- protDiaDeferredDesignWorkflow(spec)
    context <- workflow$workflow_context
    identity <- context$getIdentity()
    store <- newArtifactStore(context$getPaths(), identity$project_id)
    lineage <- result$design_result$state_manifest$active_lineage
    entry <- tail(lineage, 1L)[[1L]]
    manifest <- artifactWorkflowStateReadManifest(
        store,
        entry$manifest_relative_path
    )
    list(
        workflow = workflow,
        parity_spec = protDiaS4ParitySpecFromDigest(
            store,
            manifest,
            "PeptideQuantitativeData",
            result$exact_state_digest
        )
    )
}

runProtDiaDeferredDesignParityFork <- function(parity_spec, timeout_ms) {
    job <- parallel::mcparallel(
        tryCatch(
            runProtDiaS4ParityWorker(parity_spec),
            error = function(error) list(
                valid = FALSE,
                verifier_pid = as.integer(Sys.getpid()),
                error_class = class(error),
                error_message = conditionMessage(error),
                complete_payload_returned = FALSE
            )
        ),
        mc.set.seed = FALSE,
        silent = TRUE
    )
    deadline <- proc.time()[["elapsed"]] + as.numeric(timeout_ms) / 1000
    repeat {
        collected <- parallel::mccollect(job, wait = FALSE)
        if (!is.null(collected)) return(unname(collected)[[1L]])
        if (proc.time()[["elapsed"]] >= deadline) {
            tools::pskill(job$pid, signal = 15L)
            parallel::mccollect(job, wait = TRUE)
            protDiaArtifactAbort(
                "DIA-NN parity worker exceeded its timeout",
                "multischolar_prot_dia_s4_parity_process_failed"
            )
        }
        Sys.sleep(0.05)
    }
}

rollbackProtDiaDeferredDesignVerification <- function(spec, result) {
    workflow <- protDiaDeferredDesignWorkflow(spec)
    context <- workflow$workflow_context
    try(
        protDiaArtifactUpdateStageStatus(
            context,
            result$design_result,
            completed = FALSE
        ),
        silent = TRUE
    )
    manager <- tryCatch(newProtDiaArtifactStateManager(context), error = identity)
    if (inherits(manager, "ArtifactWorkflowState")) {
        on.exit(manager$close(), add = TRUE)
        if (isTRUE(manager$hasState("initial"))) {
            try(manager$revertToState("initial"), silent = TRUE)
        }
    }
    invisible(FALSE)
}

finalizeProtDiaDeferredDesignVerification <- function(
    spec,
    result,
    timeout_ms
) {
    prepared <- protDiaDeferredDesignParitySpec(spec, result)
    proof <- tryCatch(
        if (identical(protDiaDeferredDesignParityMode(), "fork")) {
            runProtDiaDeferredDesignParityFork(prepared$parity_spec, timeout_ms)
        } else {
            runProtDiaS4ParityProcess(prepared$parity_spec, timeout_ms)
        },
        error = identity
    )
    if (inherits(proof, "error") || !isTRUE(proof$valid) ||
        isTRUE(proof$complete_payload_returned)) {
        rollbackProtDiaDeferredDesignVerification(spec, result)
        if (inherits(proof, "error")) stop(proof)
        protDiaArtifactAbort(
            "DIA-NN deferred design parity proof is invalid",
            "multischolar_prot_dia_s4_parity_process_failed"
        )
    }
    result$design_result$hydration_verification <- proof
    workflow <- prepared$workflow
    workflow$artifact_stage_results$design <- result$design_result
    parity <- verifyProtDiaBoundedStageParity(workflow, resource_policy = NULL)
    result$readthrough_proof <- parity$proof
    result$compatibility_checkpoint <- parity$compatibility_checkpoint
    result$import_refs <- parity$import_refs
    result$design_refs <- parity$design_refs
    result$state_generation_id <- parity$state_generation_id
    result
}

runProtDiaDeferredDesignProcess <- function(spec, timeout_ms = 900000L) {
    spec <- validateProtDiaDeferredDesignSpec(spec)
    if (identical(protDiaDeferredDesignWorkerMode(), "fork")) {
        result <- runProtDiaDeferredDesignFork(spec, timeout_ms)
        status <- 0L
    } else {
        process_result <- runProtDiaDeferredDesignProcessx(spec, timeout_ms)
        result <- process_result$result
        status <- process_result$status
    }
    artifactResourceDataOnly(result, "DIA-NN design worker result")
    if (!isTRUE(result$ok) || !identical(as.integer(status), 0L) ||
        isTRUE(result$complete_payload_returned)) {
        protDiaArtifactAbort(
            paste(
                "DIA-NN design worker failed:",
                result$error_message %||% "unknown worker failure"
            ),
            "multischolar_prot_dia_design_worker_failed"
        )
    }
    finalizeProtDiaDeferredDesignVerification(spec, result, timeout_ms)
}

applyProtDiaDeferredDesignResult <- function(workflow_data, result) {
    if (!is.list(result) || !isTRUE(result$ok) ||
        isTRUE(result$complete_payload_returned)) {
        protDiaArtifactAbort(
            "DIA-NN design worker result is invalid",
            "multischolar_prot_dia_design_worker_failed"
        )
    }
    workflow_data$artifact_stage_results$design <- result$design_result
    workflow_data$data_tbl <- NULL
    workflow_data$data_cln <- NULL
    workflow_data$artifact_readthrough_refs <- list(
        import = result$import_refs,
        design = result$design_refs
    )
    workflow_data$artifact_readthrough_proof <- result$readthrough_proof
    workflow_data$artifact_compatibility_checkpoint <-
        result$compatibility_checkpoint
    previous <- workflow_data$state_manager
    workflow_data$state_manager <- newProtDiaSettledStateManager(workflow_data)
    if (inherits(previous, "ArtifactWorkflowState")) {
        try(previous$close(), silent = TRUE)
    }
    eviction <- list(
        enabled = TRUE,
        ok = TRUE,
        evicted = TRUE,
        reason = "deferred_design_worker_payload_free",
        evicted_fields = PROT_DIA_EVICT_FIELDS,
        state_generation_id = result$state_generation_id,
        state_manager_replaced = TRUE,
        complete_payload_returned = FALSE
    )
    recordProtDiaArtifactResult(workflow_data, "eviction", eviction)
    invisible(eviction)
}

runProtDiaDeferredDesignSaveFlow <- function(
    results,
    workflow_data,
    experiment_paths,
    session,
    qc_trigger = NULL
) {
    annotation_input <- protDiaDeferredAnnotationInput(workflow_data)
    if (is.null(annotation_input)) {
        protDiaArtifactAbort(
            "DIA-NN deferred design lacks bounded annotation input",
            "multischolar_missing_prot_dia_design_input"
        )
    }
    workflow_data$design_matrix <- results$design_matrix
    workflow_data$contrasts_tbl <- results$contrasts_tbl
    workflow_data$config_list <- results$config_list
    workflow_data$data_cln <- annotation_input
    on.exit(workflow_data$data_cln <- NULL, add = TRUE)
    context <- workflow_data$workflow_context
    label <- context$getStaticIdentity()$omic_label
    persist_worker <- function(...) {
        spec <- protDiaDeferredDesignSpec(
            workflow_data,
            results,
            experiment_paths,
            label
        )
        worker <- runProtDiaDeferredDesignProcess(spec)
        applyProtDiaDeferredDesignResult(workflow_data, worker)
        list(
            enabled = TRUE,
            ok = TRUE,
            committed = FALSE,
            stage_id = "design",
            deferred_worker = TRUE
        )
    }
    completeProtDesignPostCheckpoint(
        workflowData = workflow_data,
        experimentPaths = experiment_paths,
        session = session,
        qcTrigger = qc_trigger,
        successMessage = "Design matrix and contrasts saved successfully!",
        debugQcTrigger = TRUE,
        persistArtifactFn = persist_worker
    )
    TRUE
}
