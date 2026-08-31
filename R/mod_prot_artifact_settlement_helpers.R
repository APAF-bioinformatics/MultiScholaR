# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

.PROT_DIA_PARITY_WORKER_SCHEMA <- "multischolar.prot_dia_parity_worker"
.PROT_DIA_PARITY_WORKER_VERSION <- 1L
.PROT_DIA_S4_PARITY_SCHEMA <- "multischolar.prot_dia_s4_parity_worker"
.PROT_DIA_S4_PARITY_VERSION <- 1L

protDiaParityWorkerSpec <- function(
    experiment_paths,
    experiment_label,
    storage_policy,
    resource_policy = NULL,
    failure_stage = NULL
) {
    spec <- list(
        schema = .PROT_DIA_PARITY_WORKER_SCHEMA,
        schema_version = .PROT_DIA_PARITY_WORKER_VERSION,
        experiment_paths = experiment_paths,
        experiment_label = experiment_label,
        storage_policy = normalizeWorkflowStoragePolicy(storage_policy),
        resource_policy = normalizeProjectRegistryPolicy(resource_policy),
        failure_stage = failure_stage
    )
    artifactResourceDataOnly(spec, "DIA-NN parity worker specification")
    spec
}

validateProtDiaParityWorkerSpec <- function(spec) {
    valid <- is.list(spec) && identical(
        names(spec),
        c(
            "schema", "schema_version", "experiment_paths",
            "experiment_label", "storage_policy", "resource_policy",
            "failure_stage"
        )
    ) && identical(spec$schema, .PROT_DIA_PARITY_WORKER_SCHEMA) &&
        identical(
            as.integer(spec$schema_version),
            .PROT_DIA_PARITY_WORKER_VERSION
        ) && is.list(spec$experiment_paths) &&
        workflowCapabilityScalarString(spec$experiment_label) &&
        is.list(spec$storage_policy) && is.list(spec$resource_policy) &&
        (is.null(spec$failure_stage) ||
            workflowCapabilityScalarString(spec$failure_stage))
    if (!isTRUE(valid)) {
        protDiaArtifactAbort(
            "DIA-NN parity worker specification is invalid",
            "multischolar_invalid_prot_dia_parity_spec"
        )
    }
    artifactResourceDataOnly(spec, "DIA-NN parity worker specification")
    spec
}

protDiaParityWorkerFailure <- function(requested, current) {
    if (identical(requested, current)) {
        protDiaArtifactAbort(
            paste("injected DIA-NN parity failure at", current),
            "multischolar_injected_prot_dia_parity_failure"
        )
    }
    invisible(NULL)
}

runProtDiaReadthroughParityWorker <- function(spec) {
    spec <- validateProtDiaParityWorkerSpec(spec)
    protDiaParityWorkerFailure(spec$failure_stage, "before_context")
    prepared <- createProtDiaResumeContext(
        spec$experiment_paths,
        spec$experiment_label,
        spec$storage_policy
    )
    if (!isTRUE(prepared$enabled)) {
        protDiaArtifactAbort(
            "DIA-NN parity worker could not open the exact artifact project",
            "multischolar_prot_dia_parity_context_unavailable"
        )
    }
    protDiaParityWorkerFailure(spec$failure_stage, "before_hydration")
    bundle <- hydrateProtDiaResumeBundle(
        prepared$context,
        spec$resource_policy,
        retain_source_payloads = FALSE
    )
    on.exit(bundle$state_manager$close(), add = TRUE)
    protDiaParityWorkerFailure(spec$failure_stage, "after_hydration")
    proof <- newProtDiaReadthroughProof(bundle)
    checkpoint <- newProtDiaCompatibilityCheckpoint(proof)
    result <- list(
        ok = TRUE,
        worker_pid = as.integer(Sys.getpid()),
        proof = proof,
        compatibility_checkpoint = checkpoint,
        import_refs = bundle$evidence$import$refs,
        design_refs = bundle$evidence$design$refs,
        project_id = bundle$evidence$identity$project_id,
        import_run_id = bundle$evidence$import$run_id,
        design_run_id = bundle$evidence$design$run_id,
        state_generation_id = bundle$evidence$current_state$generation_id,
        complete_payload_returned = FALSE
    )
    protDiaParityWorkerFailure(spec$failure_stage, "before_result")
    artifactResourceDataOnly(result, "DIA-NN parity worker result")
    result
}

runProtDiaReadthroughParityWorkerCommand <- function(spec_path, result_path) {
    result <- tryCatch(
        runProtDiaReadthroughParityWorker(readRDS(spec_path)),
        error = \(error) list(
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

protDiaParityWorkerExpression <- function(source_tree) {
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
            "suppressMessages(attach(worker_env, name = 'multischolar-worker')); ",
            "invisible(loadNamespace('tibble')); ",
            "methods::cacheMetaData(worker_env, TRUE)"
        )
    }
    paste(
        "args <- commandArgs(trailingOnly = TRUE)",
        load_package,
        paste0(
            if (isTRUE(source_tree)) {
                "MultiScholaR:::runProtDiaReadthroughParityWorkerCommand("
            } else {
                "worker_env$runProtDiaReadthroughParityWorkerCommand("
            },
            "args[[2L]], args[[3L]])"
        ),
        sep = "; "
    )
}

runProtDiaReadthroughParityProcess <- function(spec, timeout_ms = 600000L) {
    spec <- validateProtDiaParityWorkerSpec(spec)
    spec_path <- tempfile("prot-dia-parity-spec-", fileext = ".rds")
    result_path <- tempfile("prot-dia-parity-result-", fileext = ".rds")
    on.exit(unlink(c(spec_path, result_path), force = TRUE), add = TRUE)
    saveRDS(spec, spec_path, version = 3L)
    namespace_path <- getNamespaceInfo(asNamespace("MultiScholaR"), "path")
    source_tree <- file.exists(file.path(
        namespace_path,
        "R",
        "mod_prot_artifact_settlement_helpers.R"
    ))
    process <- processx::run(
        command = file.path(R.home("bin"), "Rscript"),
        args = c(
            "--vanilla", "-e", protDiaParityWorkerExpression(source_tree),
            namespace_path, spec_path, result_path
        ),
        error_on_status = FALSE,
        timeout = as.numeric(timeout_ms) / 1000,
        echo = FALSE
    )
    if (!file.exists(result_path)) {
        protDiaArtifactAbort(
            "DIA-NN parity worker exited without a result",
            "multischolar_prot_dia_parity_process_failed"
        )
    }
    result <- readRDS(result_path)
    artifactResourceDataOnly(result, "DIA-NN parity worker result")
    if (!isTRUE(result$ok) || !identical(as.integer(process$status), 0L) ||
        isTRUE(result$complete_payload_returned)) {
        protDiaArtifactAbort(
            "DIA-NN parity worker failed",
            "multischolar_prot_dia_parity_process_failed",
            worker_error_class = result$error_class %||% character()
        )
    }
    result
}

protDiaStageParityAvailable <- function(workflow_data) {
    import <- workflow_data$artifact_stage_results$import
    design <- workflow_data$artifact_stage_results$design
    is.list(import$hydration_proofs) && is.list(import$exact_digests) &&
        is.list(design$hydration_verification)
}

validateProtDiaImportStageParity <- function(import, evidence) {
    roles <- names(evidence$import$refs)
    valid <- setequal(names(import$hydration_proofs), roles) &&
        setequal(names(import$exact_digests), roles) &&
        all(vapply(roles, function(role) {
            ref <- artifactStoreNormalizeRef(evidence$import$refs[[role]])
            proof <- import$hydration_proofs[[role]]
            identical(proof$artifact_id, ref$artifact_id) &&
                identical(
                    proof$semantic_digest,
                    ref$hash_policy$semantic$digest
                ) && identical(
                    proof$byte_digest,
                    ref$hash_policy$byte$digest
                ) && identical(
                    proof$hydration_digest,
                    import$exact_digests[[role]]
                )
        }, logical(1)))
    if (!isTRUE(valid)) {
        protDiaArtifactAbort(
            "DIA-NN import parity proof differs from committed evidence",
            "multischolar_inexact_prot_dia_stage_parity"
        )
    }
    invisible(TRUE)
}

validateProtDiaDesignStageParity <- function(design, evidence) {
    proof <- design$hydration_verification
    manifest <- artifactWorkflowStateReadManifest(
        evidence$store,
        evidence$current_state$manifest_relative_path
    )
    valid <- is.list(proof) && isTRUE(proof$valid) &&
        proof$mode %in% c("inline_exact", "fresh_process_exact") &&
        (identical(proof$mode, "inline_exact") ||
            !identical(proof$verifier_pid, as.integer(Sys.getpid()))) &&
        identical(proof$expected_digest, proof$hydrated_digest) &&
        identical(proof$generation_id, evidence$current_state$generation_id) &&
        identical(
            proof$manifest_semantic_digest,
            manifest$data$semantic_digest
        ) && !isTRUE(proof$complete_payload_returned)
    if (!isTRUE(valid)) {
        protDiaArtifactAbort(
            "DIA-NN design parity proof differs from current state evidence",
            "multischolar_inexact_prot_dia_stage_parity"
        )
    }
    invisible(TRUE)
}

verifyProtDiaBoundedStageParity <- function(workflow_data, resource_policy) {
    evidence <- collectProtDiaResumeEvidence(
        workflow_data$workflow_context,
        resource_policy,
        payload_validation = "digest"
    )
    protDiaArtifactValidateResumeContract(evidence)
    import <- workflow_data$artifact_stage_results$import
    design <- workflow_data$artifact_stage_results$design
    validateProtDiaImportStageParity(import, evidence)
    validateProtDiaDesignStageParity(design, evidence)
    bundle <- list(
        evidence = evidence,
        annotation_completed = TRUE,
        readthrough_mode = "settled"
    )
    proof <- newProtDiaReadthroughProof(bundle)
    checkpoint <- newProtDiaCompatibilityCheckpoint(proof)
    list(
        ok = TRUE,
        worker_pid = design$hydration_verification$verifier_pid,
        proof = proof,
        compatibility_checkpoint = checkpoint,
        import_refs = evidence$import$refs,
        design_refs = evidence$design$refs,
        project_id = evidence$identity$project_id,
        import_run_id = evidence$import$run_id,
        design_run_id = evidence$design$run_id,
        state_generation_id = evidence$current_state$generation_id,
        complete_payload_returned = FALSE,
        parity_reused = TRUE
    )
}

protDiaS4ParitySpec <- function(store, manifest, expected_object) {
    store <- validateArtifactStore(store)
    manifest <- artifactWorkflowStateValidateManifest(manifest)
    spec <- list(
        schema = .PROT_DIA_S4_PARITY_SCHEMA,
        schema_version = .PROT_DIA_S4_PARITY_VERSION,
        store = store,
        generation_id = manifest$generation_id,
        expected_class = if (is.null(expected_object)) {
            NULL
        } else {
            class(expected_object)[[1L]]
        },
        expected_digest = artifactExactS4HydrationDigest(expected_object)
    )
    artifactResourceDataOnly(spec, "DIA-NN S4 parity worker specification")
    spec
}

validateProtDiaS4ParitySpec <- function(spec) {
    checks <- list(
        list_value = is.list(spec),
        fields = is.list(spec) && identical(names(spec), c(
            "schema", "schema_version", "store", "generation_id",
            "expected_class", "expected_digest"
        )),
        schema = is.list(spec) && identical(
            spec$schema,
            .PROT_DIA_S4_PARITY_SCHEMA
        ),
        version = is.list(spec) && identical(
            as.integer(spec$schema_version),
            .PROT_DIA_S4_PARITY_VERSION
        ),
        store = is.list(spec) && is.list(spec$store),
        generation_id = is.list(spec) &&
            workflowCapabilityScalarString(spec$generation_id),
        expected_class = is.list(spec) &&
            (is.null(spec$expected_class) ||
                workflowCapabilityScalarString(spec$expected_class)),
        expected_digest = is.list(spec) &&
            workflowCapabilityScalarString(spec$expected_digest) &&
            grepl("^[0-9a-f]{64}$", spec$expected_digest)
    )
    if (!all(unlist(checks, use.names = FALSE))) {
        protDiaArtifactAbort(
            paste(
                "DIA-NN S4 parity worker specification is invalid:",
                paste(names(checks)[!unlist(checks, use.names = FALSE)],
                    collapse = ",")
            ),
            "multischolar_invalid_prot_dia_s4_parity_spec"
        )
    }
    artifactResourceDataOnly(spec, "DIA-NN S4 parity worker specification")
    spec
}

runProtDiaS4ParityWorker <- function(spec) {
    spec <- validateProtDiaS4ParitySpec(spec)
    manifest <- artifactWorkflowStateReadManifest(
        spec$store,
        artifactWorkflowStateManifestRelativePath(
            spec$store,
            spec$generation_id
        )
    )
    manifest_digest <- manifest$data$semantic_digest
    if (!is.character(manifest_digest) || length(manifest_digest) != 1L ||
        is.na(manifest_digest) || !grepl("^[0-9a-f]{64}$", manifest_digest)) {
        protDiaArtifactAbort(
            paste0(
                "DIA-NN S4 parity manifest digest is invalid: generation=",
                manifest$generation_id %||% "<missing>",
                ";data_fields=",
                paste(names(manifest$data), collapse = ","),
                ";digest_type=", typeof(manifest_digest),
                ";digest_length=", length(manifest_digest)
            ),
            "multischolar_invalid_prot_dia_s4_parity_spec"
        )
    }
    hydrated <- artifactWorkflowStateHydrateData(
        spec$store,
        manifest,
        hydrateDiaS4Artifact
    )
    hydrated_class <- if (is.null(hydrated)) NULL else class(hydrated)[[1L]]
    hydrated_digest <- artifactExactS4HydrationDigest(hydrated)
    valid <- identical(hydrated_class, spec$expected_class) &&
        identical(hydrated_digest, spec$expected_digest) &&
        (!isS4(hydrated) ||
            identical(methods::validObject(hydrated, test = TRUE), TRUE))
    if (!isTRUE(valid)) {
        protDiaArtifactAbort(
            paste0(
                "DIA-NN S4 parity worker found a hydration mismatch: ",
                "expected_class=", spec$expected_class %||% "<null>",
                ";hydrated_class=", hydrated_class %||% "<null>",
                ";expected_digest=", spec$expected_digest,
                ";hydrated_digest=", hydrated_digest
            ),
            "multischolar_inexact_prot_dia_s4_parity"
        )
    }
    result <- list(
        valid = TRUE,
        mode = "fresh_process_exact",
        verifier_pid = as.integer(Sys.getpid()),
        expected_digest = spec$expected_digest,
        hydrated_digest = hydrated_digest,
        manifest_semantic_digest = manifest_digest,
        generation_id = manifest$generation_id,
        complete_payload_returned = FALSE
    )
    hydrated <- NULL
    artifactResourceDataOnly(result, "DIA-NN S4 parity worker result")
    result
}

runProtDiaS4ParityWorkerCommand <- function(spec_path, result_path) {
    result <- tryCatch(
        runProtDiaS4ParityWorker(readRDS(spec_path)),
        error = \(error) list(
            valid = FALSE,
            verifier_pid = as.integer(Sys.getpid()),
            error_class = class(error),
            error_message = conditionMessage(error),
            complete_payload_returned = FALSE
        )
    )
    saveRDS(result, result_path, version = 3L)
    if (!isTRUE(result$valid)) stop(result$error_message, call. = FALSE)
    invisible(TRUE)
}

protDiaS4ParityWorkerExpression <- function(source_tree) {
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
            "suppressMessages(attach(worker_env, name = 'multischolar-worker')); ",
            "invisible(loadNamespace('tibble')); ",
            "methods::cacheMetaData(worker_env, TRUE)"
        )
    }
    paste(
        "args <- commandArgs(trailingOnly = TRUE)",
        load_package,
        paste0(
            if (isTRUE(source_tree)) {
                "MultiScholaR:::runProtDiaS4ParityWorkerCommand("
            } else {
                "worker_env$runProtDiaS4ParityWorkerCommand("
            },
            "args[[2L]], args[[3L]])"
        ),
        sep = "; "
    )
}

runProtDiaS4ParityProcess <- function(spec, timeout_ms = 600000L) {
    spec <- validateProtDiaS4ParitySpec(spec)
    spec_path <- tempfile("prot-dia-s4-parity-spec-", fileext = ".rds")
    result_path <- tempfile("prot-dia-s4-parity-result-", fileext = ".rds")
    on.exit(unlink(c(spec_path, result_path), force = TRUE), add = TRUE)
    saveRDS(spec, spec_path, version = 3L)
    namespace_path <- getNamespaceInfo(asNamespace("MultiScholaR"), "path")
    source_tree <- file.exists(file.path(
        namespace_path,
        "R",
        "mod_prot_artifact_settlement_helpers.R"
    ))
    process <- processx::run(
        command = file.path(R.home("bin"), "Rscript"),
        args = c(
            "--vanilla", "-e", protDiaS4ParityWorkerExpression(source_tree),
            namespace_path, spec_path, result_path
        ),
        error_on_status = FALSE,
        timeout = as.numeric(timeout_ms) / 1000,
        echo = FALSE
    )
    if (!file.exists(result_path)) {
        protDiaArtifactAbort(
            "DIA-NN S4 parity worker exited without a result",
            "multischolar_prot_dia_s4_parity_process_failed"
        )
    }
    result <- readRDS(result_path)
    artifactResourceDataOnly(result, "DIA-NN S4 parity worker result")
    if (!isTRUE(result$valid) || !identical(as.integer(process$status), 0L) ||
        isTRUE(result$complete_payload_returned)) {
        detail <- result$error_message %||% "unknown parity worker failure"
        protDiaArtifactAbort(
            paste("DIA-NN S4 parity worker failed:", detail),
            "multischolar_prot_dia_s4_parity_process_failed",
            worker_error_class = result$error_class %||% character(),
            worker_error_message = detail
        )
    }
    result
}

verifyProtDiaArtifactStateInWorker <- function(
    store,
    manifest,
    expected_object,
    hydrate_fn
) {
    if (!is.function(hydrate_fn)) {
        protDiaArtifactAbort(
            "DIA-NN S4 parity hydrate function is invalid",
            "multischolar_invalid_prot_dia_s4_parity_spec"
        )
    }
    if (is.null(expected_object) &&
        is.null(manifest$data$semantic_digest)) {
        proof <- list(
            valid = TRUE,
            mode = "empty_state_no_payload",
            verifier_pid = as.integer(Sys.getpid()),
            expected_digest = NULL,
            hydrated_digest = NULL,
            generation_id = manifest$generation_id,
            complete_payload_returned = FALSE
        )
        artifactResourceDataOnly(proof, "DIA-NN empty-state parity proof")
        return(proof)
    }
    runProtDiaS4ParityProcess(protDiaS4ParitySpec(
        store,
        manifest,
        expected_object
    ))
}
