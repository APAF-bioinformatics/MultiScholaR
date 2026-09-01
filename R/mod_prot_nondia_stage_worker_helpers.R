# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

.PROT_NONDIA_PENDING_STAGE_SCHEMA <-
    "multischolar.prot_nondia_pending_stage"

#' Resolve parser parameters for one non-DIA import request
#'
#' @param format Exact non-DIA input format.
#' @param parser_parameters Optional parser parameter overrides.
#'
#' @return A validated named parser parameter list.
#' @noRd
protNonDiaWorkerParserParameters <- function(
    format,
    parser_parameters = NULL
) {
    spec <- protNonDiaArtifactImportSpec(format)
    if (is.null(spec) || !format %in% c("maxquant", "fragpipe")) {
        protNonDiaArtifactAbort(
            "non-DIA import worker requires a supported LFQ format",
            "multischolar_invalid_prot_nondia_worker_format"
        )
    }
    if (is.null(parser_parameters)) parser_parameters <- list()
    if (!is.list(parser_parameters) ||
        length(setdiff(names(parser_parameters), names(spec$parser_parameters)))) {
        protNonDiaArtifactAbort(
            "non-DIA import worker parser parameters are invalid",
            "multischolar_invalid_prot_nondia_worker_parameters"
        )
    }
    parameters <- utils::modifyList(spec$parser_parameters, parser_parameters)
    valid <- all(vapply(parameters, \(value) {
        is.logical(value) && length(value) == 1L && !is.na(value)
    }, logical(1)))
    if (!isTRUE(valid)) {
        protNonDiaArtifactAbort(
            "non-DIA import worker parser parameters must be scalar flags",
            "multischolar_invalid_prot_nondia_worker_parameters"
        )
    }
    parameters
}

#' Build parser parameters from proteomics import inputs
#'
#' @param format Exact selected input format.
#' @param input Shiny input object or list.
#'
#' @return A named parser parameter list, or `NULL` outside non-DIA LFQ.
#' @noRd
protImportArtifactParserParameters <- function(format, input) {
    switch(format,
        maxquant = list(
            use_lfq = input$maxquant_use_lfq %||% TRUE,
            filter_contaminants =
                input$maxquant_filter_contaminants %||% TRUE
        ),
        fragpipe = list(
            use_maxlfq = input$fragpipe_use_maxlfq %||% TRUE
        ),
        NULL
    )
}

#' Create one non-DIA import worker specification
#'
#' @param mode Worker mode.
#' @param stage Payload-free prepared stage.
#' @param source_path Exact source file.
#' @param format Exact input format.
#' @param parser_parameters Validated parser parameters.
#' @param sanitize_names Whether run names are sanitized.
#' @param failure_stage Optional failure injection marker.
#'
#' @return A serializable worker specification.
#' @noRd
protNonDiaArtifactWorkerSpec <- function(
    mode,
    stage,
    source_path,
    format,
    parser_parameters,
    sanitize_names,
    failure_stage = NULL
) {
    list(
        mode = mode,
        stage = stage,
        source_path = normalizePath(source_path, mustWork = TRUE),
        format = format,
        parser_parameters = parser_parameters,
        sanitize_names = sanitize_names,
        failure_stage = failure_stage
    )
}

#' Validate one non-DIA import worker specification
#'
#' @param spec Candidate worker specification.
#' @param mode Required worker mode.
#'
#' @return The validated specification.
#' @noRd
validateProtNonDiaArtifactWorkerSpec <- function(spec, mode) {
    valid <- is.list(spec) && identical(spec$mode, mode) &&
        mode %in% c("writer", "verifier") &&
        is.list(spec$stage) && identical(spec$stage$stage_id, "import") &&
        workflowCapabilityScalarString(spec$source_path) &&
        file.exists(spec$source_path) &&
        spec$format %in% c("maxquant", "fragpipe") &&
        is.logical(spec$sanitize_names) && length(spec$sanitize_names) == 1L &&
        !is.na(spec$sanitize_names)
    if (!isTRUE(valid)) {
        protNonDiaArtifactAbort(
            "non-DIA import worker specification is invalid",
            "multischolar_invalid_prot_nondia_worker_spec"
        )
    }
    spec$parser_parameters <- protNonDiaWorkerParserParameters(
        spec$format,
        spec$parser_parameters
    )
    spec
}

#' Inject one non-DIA worker failure
#'
#' @param requested_stage Requested failure marker.
#' @param current_stage Current worker marker.
#'
#' @return `TRUE`, invisibly, unless the failure is requested.
#' @noRd
protNonDiaArtifactWorkerFailure <- function(requested_stage, current_stage) {
    if (identical(requested_stage, current_stage)) {
        protNonDiaArtifactAbort(
            paste("injected non-DIA worker failure at", current_stage),
            "multischolar_injected_prot_nondia_worker_failure"
        )
    }
    invisible(TRUE)
}

#' Run the reviewed parser for one non-DIA LFQ worker
#'
#' @param spec Validated worker specification.
#'
#' @return The established public importer result.
#' @noRd
protNonDiaArtifactWorkerImport <- function(spec) {
    imported <- switch(spec$format,
        maxquant = do.call(
            importMaxQuantData,
            c(list(filepath = spec$source_path), spec$parser_parameters)
        ),
        fragpipe = do.call(
            importFragPipeData,
            c(list(filepath = spec$source_path), spec$parser_parameters)
        )
    )
    if (isTRUE(spec$sanitize_names)) {
        imported <- sanitizeProtDiaArtifactImport(imported)
    }
    imported
}

#' Run one non-DIA import writer worker
#'
#' @param spec Writer specification.
#'
#' @return A payload-free writer result.
#' @noRd
runProtNonDiaArtifactWriterWorker <- function(spec) {
    spec <- validateProtNonDiaArtifactWorkerSpec(spec, "writer")
    protNonDiaArtifactWorkerFailure(spec$failure_stage, "before_import")
    output_path <- tempfile("prot-nondia-streaming-", fileext = ".parquet")
    on.exit(unlink(output_path, force = TRUE), add = TRUE)
    streamed <- writeProtNonDiaStreamingParquet(
        spec$source_path,
        output_path,
        spec$format,
        spec$parser_parameters
    )
    encoded <- encodeArtifactStreamingParquet(
        output_path,
        owner = paste(
            spec$stage$descriptor$descriptor_id,
            "import",
            "canonical_data",
            sep = "."
        )
    )
    protNonDiaArtifactWorkerFailure(spec$failure_stage, "after_import")
    ref <- artifactStorePublishStreamingParquet(
        spec$stage$store,
        encoded,
        logical_key = list(
            project_id = spec$stage$identity$project_id,
            omic_type = spec$stage$identity$omic_type,
            workflow_slug = spec$stage$identity$workflow_slug,
            stage_id = "import",
            state_role = "canonical_data",
            generation_id = spec$stage$generation_id
        ),
        provenance_ids = c(spec$stage$run_id, spec$stage$action_id),
        verification = "deferred_exact",
    )
    protNonDiaArtifactWorkerFailure(spec$failure_stage, "after_write")
    if (identical(spec$failure_stage, "hard_exit_after_write")) {
        quit(save = "no", status = 70L, runLast = FALSE)
    }
    result <- list(
        ok = TRUE,
        mode = "writer",
        worker_pid = as.integer(Sys.getpid()),
        refs = list(canonical_data = unclass(ref)),
        exact_digests = list(canonical_data = NULL),
        data_type = "protein",
        column_mapping = streamed$column_mapping,
        source = protNonDiaArtifactSourceMetadata(
            spec$source_path,
            spec$format
        ),
        bounded_streaming = TRUE
    )
    protNonDiaArtifactWorkerFailure(spec$failure_stage, "before_result")
    artifactResourceDataOnly(result, "non-DIA import writer result")
    result
}

#' Run one independent non-DIA import verifier worker
#'
#' @param spec Verifier specification.
#'
#' @return A payload-free exact parity result.
#' @noRd
runProtNonDiaArtifactVerifierWorker <- function(spec) {
    spec <- validateProtNonDiaArtifactWorkerSpec(spec, "verifier")
    protNonDiaArtifactWorkerFailure(spec$failure_stage, "before_verify")
    proofs <- lapply(spec$stage$refs, \(ref) {
        artifactStoreVerifyStreamingRef(
            spec$stage$store,
            ref,
            verify_semantic = FALSE
        )
    })
    names(proofs) <- names(spec$stage$refs)
    oracle_digest <- spec$stage$refs$canonical_data$hash_policy$semantic$digest
    protNonDiaArtifactWorkerFailure(spec$failure_stage, "after_verify")
    result <- list(
        ok = TRUE,
        mode = "verifier",
        worker_pid = as.integer(Sys.getpid()),
        proofs = proofs,
        oracle_digest = oracle_digest,
        bounded_streaming = TRUE
    )
    artifactResourceDataOnly(result, "non-DIA import verifier result")
    result
}

#' Execute one non-DIA stage worker command
#'
#' @param spec_path Serialized worker specification path.
#' @param result_path Serialized result path.
#' @param mode Worker mode.
#'
#' @return `TRUE`, invisibly.
#' @noRd
runProtNonDiaArtifactStageWorkerCommand <- function(
    spec_path,
    result_path,
    mode
) {
    result <- tryCatch(
        {
            spec <- readRDS(spec_path)
            switch(mode,
                writer = runProtNonDiaArtifactWriterWorker(spec),
                verifier = runProtNonDiaArtifactVerifierWorker(spec),
                protNonDiaArtifactAbort(
                    "non-DIA artifact worker mode is unsupported",
                    "multischolar_invalid_prot_nondia_worker_mode"
                )
            )
        },
        error = \(error) list(
            ok = FALSE,
            mode = mode,
            worker_pid = as.integer(Sys.getpid()),
            error_class = class(error),
            error_message = conditionMessage(error)
        )
    )
    saveRDS(result, result_path, version = 3L)
    if (!isTRUE(result$ok)) stop(result$error_message, call. = FALSE)
    invisible(TRUE)
}

#' Build the standalone non-DIA worker expression
#'
#' @param source_tree Whether the package is loaded from a source tree.
#'
#' @return An R expression string.
#' @noRd
protNonDiaArtifactWorkerExpression <- function(source_tree) {
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
    call <- if (isTRUE(source_tree)) {
        "MultiScholaR:::runProtNonDiaArtifactStageWorkerCommand("
    } else {
        "worker_env$runProtNonDiaArtifactStageWorkerCommand("
    }
    paste(
        "args <- commandArgs(trailingOnly = TRUE)",
        load_package,
        paste0(call, "args[[2L]], args[[3L]], args[[4L]])"),
        sep = "; "
    )
}

#' Run one non-DIA import worker in a fresh process
#'
#' @param spec Worker specification.
#' @param timeout_ms Worker timeout in milliseconds.
#'
#' @return A payload-free worker result.
#' @noRd
runProtNonDiaArtifactStageProcessx <- function(spec, timeout_ms = 600000L) {
    if (!requireNamespace("processx", quietly = TRUE)) {
        protNonDiaArtifactAbort(
            "non-DIA artifact workers require package 'processx'",
            "multischolar_missing_prot_nondia_worker_dependency"
        )
    }
    spec <- validateProtNonDiaArtifactWorkerSpec(spec, spec$mode)
    spec_path <- tempfile("prot-nondia-worker-spec-", fileext = ".rds")
    result_path <- tempfile("prot-nondia-worker-result-", fileext = ".rds")
    on.exit(unlink(c(spec_path, result_path), force = TRUE), add = TRUE)
    saveRDS(spec, spec_path, version = 3L)
    namespace_path <- getNamespaceInfo(asNamespace("MultiScholaR"), "path")
    source_tree <- file.exists(file.path(
        namespace_path,
        "R",
        "mod_prot_nondia_stage_worker_helpers.R"
    ))
    started <- proc.time()[["elapsed"]]
    process <- processx::run(
        command = file.path(R.home("bin"), "Rscript"),
        args = c(
            "--vanilla", "-e",
            protNonDiaArtifactWorkerExpression(source_tree),
            namespace_path, spec_path, result_path, spec$mode
        ),
        error_on_status = FALSE,
        timeout = as.numeric(timeout_ms) / 1000,
        echo = FALSE
    )
    elapsed_seconds <- proc.time()[["elapsed"]] - started
    if (!file.exists(result_path)) {
        protNonDiaArtifactAbort(
            sprintf("non-DIA artifact %s worker returned no result", spec$mode),
            "multischolar_prot_nondia_worker_process_failed",
            worker_status = as.integer(process$status)
        )
    }
    result <- readRDS(result_path)
    result$process_elapsed_seconds <- unname(elapsed_seconds)
    artifactResourceDataOnly(result, "non-DIA import worker result")
    if (!isTRUE(result$ok) || !identical(as.integer(process$status), 0L)) {
        protNonDiaArtifactAbort(
            sprintf("non-DIA artifact %s worker failed", spec$mode),
            "multischolar_prot_nondia_worker_process_failed",
            worker_status = as.integer(process$status),
            worker_error_class = result$error_class %||% character(),
            worker_error_message = result$error_message %||% NULL
        )
    }
    result
}

#' Run one non-DIA import worker in a fork child
#'
#' @param spec Worker specification.
#' @param timeout_ms Worker timeout in milliseconds.
#'
#' @return A payload-free worker result.
#' @noRd
runProtNonDiaArtifactStageFork <- function(spec, timeout_ms = 600000L) {
    spec <- validateProtNonDiaArtifactWorkerSpec(spec, spec$mode)
    job <- parallel::mcparallel(
        tryCatch(
            switch(spec$mode,
                writer = runProtNonDiaArtifactWriterWorker(spec),
                verifier = runProtNonDiaArtifactVerifierWorker(spec)
            ),
            error = \(error) list(
                ok = FALSE,
                mode = spec$mode,
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
                paste("non-DIA", spec$mode, "worker exceeded its timeout"),
                "multischolar_prot_nondia_worker_process_failed"
            )
        }
        Sys.sleep(0.05)
    }
    if (!is.list(result) || !isTRUE(result$ok)) {
        error_class <- if (is.list(result)) {
            result$error_class %||% character()
        } else {
            class(result)
        }
        protNonDiaArtifactAbort(
            paste("non-DIA artifact", spec$mode, "worker failed"),
            unique(c(
                error_class,
                "multischolar_prot_nondia_worker_process_failed"
            )),
            worker_error_message = if (is.list(result)) {
                result$error_message %||% NULL
            } else {
                as.character(result)
            }
        )
    }
    artifactResourceDataOnly(result, "non-DIA fork worker result")
    result
}

#' Resolve the non-DIA import worker process model
#'
#' @return Either `"fork"` on Linux or `"process"` elsewhere.
#' @noRd
protNonDiaArtifactWorkerMode <- function() {
    requested <- getOption(
        "multischolar.prot_nondia.import_worker_mode",
        if (identical(Sys.info()[["sysname"]], "Linux")) "fork" else "process"
    )
    match.arg(requested, c("fork", "process"))
}

#' Run one non-DIA import worker with the platform process model
#'
#' @param spec Worker specification.
#' @param timeout_ms Worker timeout in milliseconds.
#'
#' @return A payload-free worker result.
#' @noRd
runProtNonDiaArtifactStageProcess <- function(spec, timeout_ms = 600000L) {
    if (identical(protNonDiaArtifactWorkerMode(), "fork")) {
        runProtNonDiaArtifactStageFork(spec, timeout_ms)
    } else {
        runProtNonDiaArtifactStageProcessx(spec, timeout_ms)
    }
}

#' Prepare a payload-free non-DIA worker stage
#'
#' @param prepared Prepared exact-tuple context.
#' @param resource_policy Optional registry resource policy.
#'
#' @return A serializable stage resource bundle.
#' @noRd
protNonDiaArtifactPreparedWorkerStage <- function(
    prepared,
    resource_policy = NULL
) {
    stage <- artifactStageResources(
        prepared$context,
        prepared$descriptor,
        "import",
        resource_policy,
        protNonDiaArtifactAbort
    )
    closeProjectRegistry(stage$session)
    stage$session <- NULL
    stage$resource_policy <- normalizeProjectRegistryPolicy(resource_policy)
    artifactResourceDataOnly(stage, "non-DIA pending import stage")
    stage
}

#' Build and validate a pending non-DIA import stage
#'
#' @param stage Prepared stage resources.
#' @param writer Writer result.
#' @param verifier Verifier result.
#' @param parser_parameters Parser parameters.
#' @param sanitize_names Whether run names were sanitized.
#'
#' @return A payload-free pending stage.
#' @noRd
newProtNonDiaArtifactPendingStage <- function(
    stage,
    writer,
    verifier,
    parser_parameters,
    sanitize_names
) {
    stage$refs <- writer$refs
    stage$exact_digests <- list(canonical_data = verifier$oracle_digest)
    pending <- list(
        schema = .PROT_NONDIA_PENDING_STAGE_SCHEMA,
        schema_version = 1L,
        stage = stage,
        parameters = list(
            capability_id = stage$descriptor$descriptor_id,
            column_mapping = writer$column_mapping,
            column_mapping_serialized =
                artifactWorkflowStateSerializeMetadata(
                    writer$column_mapping,
                    "non-DIA import column mapping"
                ),
            readthrough_contract_version = 1L,
            workflow_type = stage$descriptor$identity$workflow_type,
            input_format = stage$descriptor$identity$input_format,
            data_level = "protein",
            parser_parameters = parser_parameters,
            sanitize_names = isTRUE(sanitize_names)
        ),
        source = writer$source,
        proofs = verifier$proofs,
        process_evidence = list(
            writer_pid = writer$worker_pid,
            verifier_pid = verifier$worker_pid,
            distinct_workers = !identical(
                writer$worker_pid,
                verifier$worker_pid
            ),
            complete_payload_returned = FALSE,
            oracle_digest = verifier$oracle_digest
        )
    )
    validateProtNonDiaArtifactPendingStage(pending)
}

#' Validate one pending non-DIA import stage
#'
#' @param pending Candidate pending stage.
#'
#' @return The validated payload-free stage.
#' @noRd
validateProtNonDiaArtifactPendingStage <- function(pending) {
    valid <- is.list(pending) &&
        identical(pending$schema, .PROT_NONDIA_PENDING_STAGE_SCHEMA) &&
        identical(as.integer(pending$schema_version), 1L) &&
        is.list(pending$stage) && identical(pending$stage$stage_id, "import") &&
        identical(names(pending$stage$refs), "canonical_data") &&
        identical(names(pending$proofs), "canonical_data") &&
        isTRUE(pending$process_evidence$distinct_workers) &&
        identical(pending$process_evidence$complete_payload_returned, FALSE) &&
        identical(
            pending$process_evidence$oracle_digest,
            pending$stage$exact_digests$canonical_data
        )
    if (!isTRUE(valid)) {
        protNonDiaArtifactAbort(
            "pending non-DIA import stage is malformed or unverified",
            "multischolar_invalid_prot_nondia_pending_stage"
        )
    }
    artifactResourceDataOnly(pending, "non-DIA pending import stage")
    pending
}

#' Stage one MaxQuant or FragPipe import outside the parent process
#'
#' @param workflow_data Mutable proteomics workflow state.
#' @param source_path Exact source file.
#' @param format Exact input format.
#' @param parser_parameters Parser parameters.
#' @param sanitize_names Whether run names are sanitized.
#' @param resource_policy Optional registry resource policy.
#' @param timeout_ms Worker timeout in milliseconds.
#' @param writer_failure_stage Optional writer failure marker.
#' @param verifier_failure_stage Optional verifier failure marker.
#'
#' @return A staged import result and payload-free pending generation.
#' @noRd
stageProtNonDiaImportArtifacts <- function(
    workflow_data,
    source_path,
    format,
    parser_parameters = NULL,
    sanitize_names = FALSE,
    resource_policy = NULL,
    timeout_ms = 600000L,
    writer_failure_stage = NULL,
    verifier_failure_stage = NULL
) {
    if (!format %in% c("maxquant", "fragpipe")) {
        return(list(
            enabled = FALSE,
            attempted = FALSE,
            ok = TRUE,
            reason = "not_supported_nondia_lfq_worker"
        ))
    }
    parameters <- protNonDiaWorkerParserParameters(format, parser_parameters)
    prepared <- prepareProtNonDiaArtifactContext(
        workflow_data,
        format = format,
        data_type = "protein"
    )
    if (!isTRUE(prepared$enabled)) {
        return(list(
            enabled = FALSE,
            attempted = FALSE,
            ok = TRUE,
            reason = prepared$reason
        ))
    }
    stage <- protNonDiaArtifactPreparedWorkerStage(prepared, resource_policy)
    caller_owns_stage <- FALSE
    on.exit({
        if (!caller_owns_stage) {
            try(discardProtDiaArtifactWorkerStage(stage), silent = TRUE)
        }
    }, add = TRUE)
    writer <- runProtNonDiaArtifactStageProcess(
        protNonDiaArtifactWorkerSpec(
            "writer", stage, source_path, format, parameters,
            sanitize_names, writer_failure_stage
        ),
        timeout_ms
    )
    verifier_stage <- stage
    verifier_stage$refs <- writer$refs
    verifier_stage$exact_digests <- writer$exact_digests
    verifier <- runProtNonDiaArtifactStageProcess(
        protNonDiaArtifactWorkerSpec(
            "verifier", verifier_stage, source_path, format, parameters,
            sanitize_names, verifier_failure_stage
        ),
        timeout_ms
    )
    pending <- newProtNonDiaArtifactPendingStage(
        stage,
        writer,
        verifier,
        parameters,
        sanitize_names
    )
    imported <- protNonDiaArtifactWorkerImport(protNonDiaArtifactWorkerSpec(
        "writer",
        stage,
        source_path,
        format,
        parameters,
        sanitize_names
    ))
    ref <- pending$stage$refs$canonical_data
    shape_valid <- identical(
        as.numeric(ref$shape$rows),
        as.numeric(nrow(imported$data))
    ) &&
        identical(as.integer(ref$shape$columns), as.integer(ncol(imported$data)))
    if (!isTRUE(shape_valid)) {
        protNonDiaArtifactAbort(
            "streamed non-DIA shape differs from the reviewed importer",
            "multischolar_inexact_prot_nondia_worker_hydration"
        )
    }
    if (!identical(imported$column_mapping, writer$column_mapping)) {
        protNonDiaArtifactAbort(
            "streamed non-DIA mapping differs from the reviewed importer",
            "multischolar_inexact_prot_nondia_import_metadata"
        )
    }
    pending$process_evidence$parser_parent_pid <- as.integer(Sys.getpid())
    list_result <- list(
        enabled = TRUE,
        attempted = TRUE,
        ok = TRUE,
        result = imported,
        pending_stage = pending,
        process_evidence = pending$process_evidence
    )
    caller_owns_stage <- TRUE
    list_result
}

#' Safely stage one non-DIA LFQ artifact import
#'
#' @inheritParams stageProtNonDiaImportArtifacts
#' @param log_warn Warning logger used for memory fallback.
#'
#' @return A staged result or explicit memory fallback evidence.
#' @noRd
stageProtNonDiaImportArtifactsSafely <- function(
    workflow_data,
    source_path,
    format,
    parser_parameters = NULL,
    sanitize_names = FALSE,
    resource_policy = NULL,
    timeout_ms = 600000L,
    writer_failure_stage = NULL,
    verifier_failure_stage = NULL,
    log_warn = logger::log_warn
) {
    result <- tryCatch(
        stageProtNonDiaImportArtifacts(
            workflow_data,
            source_path,
            format,
            parser_parameters,
            sanitize_names,
            resource_policy,
            timeout_ms,
            writer_failure_stage,
            verifier_failure_stage
        ),
        error = \(error) {
            log_warn(paste(
                "non-DIA artifact import worker failed; using memory:",
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
    audit_result <- result
    audit_result$result <- NULL
    recordArtifactStageResult(workflow_data, "import_worker", audit_result)
    result
}

#' Dispatch proteomics import staging by exact tuple
#'
#' @param workflow_data Mutable proteomics workflow state.
#' @param source_path Exact source file.
#' @param format Exact selected format.
#' @param use_precursor_norm DIA-NN precursor normalization flag.
#' @param sanitize_names Whether run names are sanitized.
#' @param parser_parameters Optional non-DIA parser parameters.
#' @param log_warn Warning logger used for memory fallback.
#' @param ... Additional worker controls used by focused tests.
#'
#' @return One exact-tuple staged import result.
#' @noRd
stageProtImportArtifactsSafely <- function(
    workflow_data,
    source_path,
    format,
    use_precursor_norm = TRUE,
    sanitize_names = FALSE,
    parser_parameters = NULL,
    log_warn = logger::log_warn,
    ...
) {
    if (identical(format, "diann")) {
        return(stageProtDiaImportArtifactsSafely(
            workflow_data,
            source_path,
            format,
            use_precursor_norm,
            sanitize_names,
            log_warn = log_warn,
            ...
        ))
    }
    stageProtNonDiaImportArtifactsSafely(
        workflow_data,
        source_path,
        format,
        parser_parameters,
        sanitize_names,
        log_warn = log_warn,
        ...
    )
}

#' Publish one verified non-DIA worker generation
#'
#' @param workflow_data Mutable proteomics workflow state.
#' @param pending_stage Verified pending import stage.
#'
#' @return A committed import stage record.
#' @noRd
publishProtNonDiaArtifactWorkerStage <- function(
    workflow_data,
    pending_stage
) {
    pending <- validateProtNonDiaArtifactPendingStage(pending_stage)
    context <- workflow_data$workflow_context
    if (!inherits(context, "WorkflowContext") || !context$isBound() ||
        !identical(context$getIdentity(), pending$stage$identity)) {
        protNonDiaArtifactAbort(
            "pending non-DIA stage belongs to another workflow context",
            "multischolar_invalid_prot_nondia_pending_stage_context"
        )
    }
    registry <- projectRegistryForContext(
        context,
        pending$stage$resource_policy
    )
    session <- initializeProjectRegistry(registry)
    on.exit(closeProjectRegistry(session), add = TRUE)
    artifactWorkflowStateEnsureWorkflow(
        session,
        pending$stage$registry_identity
    )
    stage <- pending$stage
    stage$session <- session
    artifactStageRegister(
        stage,
        stage$refs,
        pending$parameters,
        pending$source,
        deferred_commit = FALSE
    )
    list(
        enabled = TRUE,
        ok = TRUE,
        stage_id = "import",
        run_id = stage$run_id,
        action_id = stage$action_id,
        generation_id = stage$generation_id,
        refs = stage$refs,
        committed = TRUE,
        process_evidence = pending$process_evidence,
        hydration_proofs = pending$proofs,
        exact_digests = stage$exact_digests
    )
}
