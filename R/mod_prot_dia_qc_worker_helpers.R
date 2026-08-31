# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

.PROT_DIA_QC_WORKER_SCHEMA <- "multischolar.prot_dia_qc_worker"
.PROT_DIA_QC_WORKER_VERSION <- 2L

protDiaQcWorkerFunctions <- function() {
    c(
        qvalue_filter = "runPeptideQvalueApplyStep",
        intensity_filter = "runPeptideIntensityApplyStep",
        sample_filter = "runPeptideSampleApplyStep",
        protein_evidence_filter = "runProteinPeptideApplyStep",
        replicate_filter = "runPeptideReplicateApplyStep",
        precursor_rollup = "runPeptideRollupApplyStep",
        imputation = "runPeptideImputationStep",
        protein_iq_rollup = "runProteinIqRollupApplyStep",
        protein_limpa_rollup = "runProteinLimpaRollupApplyStep",
        protein_accession_cleanup = "runProteinAccessionCleanupStep",
        protein_intensity_filter = "runProteinIntensityFilterApplyStep",
        protein_duplicate_aggregation = "runProteinDuplicateRemovalStep",
        protein_replicate_filter = "runProteinReplicateFilterApplyStep"
    )
}

protDiaQcWorkerActions <- function() {
    c(names(protDiaQcWorkerFunctions()), "revert")
}

protDiaQcPlotStepNames <- function() {
    c(
        qvalue_filter = "2_qval_filtered",
        protein_evidence_filter = "3_protein_peptide_filtered",
        precursor_rollup = "4_precursor_rollup",
        intensity_filter = "5_intensity_filtered",
        sample_filter = "6_sample_filtered",
        replicate_filter = "7_replicate_filtered",
        imputation = "8_imputed",
        protein_iq_rollup = "9_protein_s4_created",
        protein_limpa_rollup = "9_protein_s4_created",
        protein_accession_cleanup = "10_protein_accession_cleaned",
        protein_intensity_filter = "11_protein_intensity_filtered",
        protein_duplicate_aggregation = "12_duplicates_removed",
        protein_replicate_filter = "13_protein_replicate_filtered"
    )
}

protDiaQcStatePlotStepNames <- function() {
    c(
        raw_data_s4 = "1_initial_data",
        qvalue_filtered = "2_qval_filtered",
        protein_peptide_filtered = "3_protein_peptide_filtered",
        precursor_rollup = "4_precursor_rollup",
        intensity_filtered = "5_intensity_filtered",
        sample_filtered = "6_sample_filtered",
        replicate_filtered = "7_replicate_filtered",
        imputed = "8_imputed",
        protein_s4_created = "9_protein_s4_created",
        protein_accession_cleaned = "10_protein_accession_cleaned",
        protein_intensity_filtered = "11_protein_intensity_filtered",
        duplicates_removed = "12_duplicates_removed",
        protein_replicate_filtered = "13_protein_replicate_filtered"
    )
}

protDiaQcWorkerEligible <- function(workflow_data) {
    context <- tryCatch(workflow_data$workflow_context, error = function(...) NULL)
    inherits(workflow_data$state_manager, "ArtifactWorkflowState") &&
        inherits(context, "WorkflowContext") && context$isBound() &&
        identical(context$getStorageDecision()$effective_rollout, "evict") &&
        !isTRUE(workflow_data$artifact_worker_active) &&
        protDiaArtifactCoordinatorOwned(workflow_data)
}

protDiaWorkerContextSpec <- function(workflow_data) {
    context <- workflow_data$workflow_context
    list(
        project_root = context$getProjectRoot(),
        static_identity = context$getStaticIdentity(),
        storage_policy = context$getStoragePolicy(),
        identity = context$getIdentity(),
        resolution = unclass(context$getStorageDecision()),
        paths = context$getPaths()
    )
}

protDiaQcWorkerSequenceRef <- function(workflow_data) {
    ref <- workflow_data$artifact_readthrough_refs$design$sequences
    if (is.list(ref)) return(ref)
    ref <- workflow_data$artifact_stage_results$design$refs$sequences
    if (is.list(ref)) ref else NULL
}

protDiaQcWorkerSessionInputs <- function(workflow_data) {
    context <- workflow_data$protein_qc_context
    sequence_ref <- protDiaQcWorkerSequenceRef(workflow_data)
    list(
        experiment_paths = context$experiment_paths,
        omic_type = context$omic_type,
        experiment_label = context$experiment_label,
        fasta_metadata = workflow_data$fasta_metadata,
        sequence_ref = sequence_ref,
        sequence_available = !is.null(workflow_data$aa_seq_tbl_final) ||
            is.list(sequence_ref)
    )
}

protDiaQcWorkerSpec <- function(workflow_data, action_id, arguments) {
    functions <- protDiaQcWorkerFunctions()
    if (!action_id %in% protDiaQcWorkerActions() || !is.list(arguments)) {
        protDiaArtifactAbort(
            "DIA-NN QC worker action is invalid",
            "multischolar_invalid_prot_dia_qc_worker_spec"
        )
    }
    spec <- list(
        schema = .PROT_DIA_QC_WORKER_SCHEMA,
        schema_version = .PROT_DIA_QC_WORKER_VERSION,
        context = protDiaWorkerContextSpec(workflow_data),
        action_id = action_id,
        arguments = arguments,
        qc_params = workflow_data$qc_params,
        session_inputs = protDiaQcWorkerSessionInputs(workflow_data)
    )
    artifactResourceDataOnly(spec, "DIA-NN QC worker specification")
    spec
}

validateProtDiaQcWorkerSpec <- function(spec) {
    required <- c(
        "schema", "schema_version", "context", "action_id", "arguments",
        "qc_params", "session_inputs"
    )
    valid <- is.list(spec) && identical(names(spec), required) &&
        identical(spec$schema, .PROT_DIA_QC_WORKER_SCHEMA) &&
        identical(as.integer(spec$schema_version), .PROT_DIA_QC_WORKER_VERSION) &&
        is.list(spec$context) && is.list(spec$arguments) &&
        is.list(spec$session_inputs) &&
        spec$action_id %in% protDiaQcWorkerActions()
    if (!isTRUE(valid)) {
        protDiaArtifactAbort(
            "DIA-NN QC worker specification is malformed",
            "multischolar_invalid_prot_dia_qc_worker_spec"
        )
    }
    artifactResourceDataOnly(spec, "DIA-NN QC worker specification")
    spec
}

protDiaQcWorkerSequenceData <- function(context, session_inputs) {
    ref <- session_inputs$sequence_ref
    if (!is.list(ref)) {
        if (isTRUE(session_inputs$sequence_available)) {
            protDiaArtifactAbort(
                "DIA-NN QC worker lacks the declared sequence artifact",
                "multischolar_missing_prot_dia_qc_worker_input"
            )
        }
        return(NULL)
    }
    identity <- context$getIdentity()
    store <- newArtifactStore(context$getPaths(), identity$project_id)
    protDiaArtifactReadTable(store, ref)
}

protDiaQcWorkerWorkflow <- function(spec) {
    value <- spec$context
    context <- WorkflowContext$new(
        value$project_root,
        value$static_identity,
        value$storage_policy
    )
    context$bind(value$identity, value$resolution, value$paths)
    workflow <- new.env(parent = emptyenv())
    workflow$workflow_context <- context
    workflow$state_manager <- newProtDiaArtifactStateManager(context)
    workflow$artifact_worker_active <- TRUE
    workflow$artifact_worker_force_materialization <- TRUE
    workflow$artifact_stage_results <- list()
    workflow$data_tbl <- NULL
    workflow$data_cln <- NULL
    workflow$data_format <- "diann"
    workflow$data_type <- "peptide"
    workflow$qc_params <- spec$qc_params
    workflow$config_list <- workflow$state_manager$getStateConfig()
    workflow$fasta_metadata <- spec$session_inputs$fasta_metadata
    workflow$aa_seq_tbl_final <- protDiaQcWorkerSequenceData(
        context,
        spec$session_inputs
    )
    workflow$protein_qc_context <- list(
        experiment_paths = spec$session_inputs$experiment_paths,
        omic_type = spec$session_inputs$omic_type,
        experiment_label = spec$session_inputs$experiment_label,
        progress_env = new.env(parent = emptyenv())
    )
    assign("config_list", workflow$config_list, envir = .GlobalEnv)
    workflow
}

protDiaQcResultObject <- function(result) {
    matches <- vapply(result, function(value) {
        isS4(value) && (
            methods::is(value, "PeptideQuantitativeData") ||
                methods::is(value, "ProteinQuantitativeData")
        )
    }, logical(1))
    if (sum(matches) != 1L) {
        protDiaArtifactAbort(
            "DIA-NN QC worker result lacks one omics S4 state",
            "multischolar_prot_dia_qc_worker_failed"
        )
    }
    result[[which(matches)]]
}

protDiaQcPlotPng <- function(object, step_name) {
    path <- tempfile("prot-dia-qc-", fileext = ".png")
    on.exit(unlink(path, force = TRUE), add = TRUE)
    grDevices::png(path, width = 1200L, height = 900L, res = 120)
    device_open <- TRUE
    on.exit({
        if (device_open) grDevices::dev.off()
    }, add = TRUE)
    data <- protDiaProteinQcTable(object)
    if (!is.data.frame(data)) {
        protDiaArtifactAbort(
            "DIA-NN QC worker cannot render its scientific state",
            "multischolar_prot_dia_qc_worker_failed"
        )
    }
    grid <- updateProteinFiltering(
        data,
        step_name,
        return_grid = TRUE,
        formats = character(),
        progress_env = new.env(parent = emptyenv())
    )
    grid::grid.newpage()
    grid::grid.draw(grid)
    grDevices::dev.off()
    device_open <- FALSE
    readBin(path, what = "raw", n = file.info(path)$size)
}

runProtDiaQcWorker <- function(spec) {
    spec <- validateProtDiaQcWorkerSpec(spec)
    workflow <- protDiaQcWorkerWorkflow(spec)
    on.exit(workflow$state_manager$close(), add = TRUE)
    if (identical(spec$action_id, "revert")) {
        state_name <- spec$arguments$state_name
        if (!workflowStateScalarString(state_name)) {
            protDiaArtifactAbort(
                "DIA-NN QC revert target is invalid",
                "multischolar_invalid_prot_dia_qc_worker_spec"
            )
        }
        object <- workflow$state_manager$revertToState(state_name)
        result <- list(
            previousState = state_name,
            resultText = paste("Reverted to previous state:", state_name)
        )
    } else {
        function_name <- protDiaQcWorkerFunctions()[[spec$action_id]]
        operation <- get(function_name, envir = environment(runProtDiaQcWorker))
        workflow_argument <- if ("workflowData" %in% names(formals(operation))) {
            "workflowData"
        } else {
            "workflow_data"
        }
        owner <- stats::setNames(list(workflow), workflow_argument)
        result <- do.call(
            operation,
            c(owner, spec$arguments)
        )
        object <- protDiaQcResultObject(result)
        result[vapply(result, isS4, logical(1))] <- NULL
    }
    plot_step <- protDiaQcStatePlotStepNames()[[
        workflow$state_manager$getCurrentStateName()
    ]]
    if (is.null(plot_step)) {
        plot_step <- protDiaQcPlotStepNames()[[spec$action_id]]
    }
    if (is.null(plot_step)) {
        protDiaArtifactAbort(
            "DIA-NN QC worker lacks a plot contract for its current state",
            "multischolar_prot_dia_qc_worker_failed"
        )
    }
    plot_png <- protDiaQcPlotPng(object, plot_step)
    output <- list(
        ok = TRUE,
        worker_pid = as.integer(Sys.getpid()),
        action_id = spec$action_id,
        result = result,
        plot_png = plot_png,
        qc_params = workflow$qc_params,
        state_export = workflow$state_manager$exportState(),
        current_state = workflow$state_manager$getCurrentStateName(),
        complete_payload_returned = FALSE
    )
    artifactResourceDataOnly(output, "DIA-NN QC worker result")
    output
}

runProtDiaQcRevert <- function(workflow_data, state_name) {
    worker <- runProtDiaQcProcess(protDiaQcWorkerSpec(
        workflow_data,
        "revert",
        list(state_name = state_name)
    ))
    result <- applyProtDiaQcWorkerResult(workflow_data, worker)
    result$plot_png <- worker$plot_png
    result$revertedS4 <- NULL
    result
}

runProtDiaQcWorkerCommand <- function(spec_path, result_path) {
    result <- tryCatch(
        runProtDiaQcWorker(readRDS(spec_path)),
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

protDiaQcWorkerExpression <- function(source_tree) {
    load <- if (isTRUE(source_tree)) {
        "pkgload::load_all(args[[1L]], quiet = TRUE, export_all = FALSE)"
    } else {
        "worker_env <- loadNamespace('MultiScholaR', lib.loc = dirname(args[[1L]]))"
    }
    call <- if (isTRUE(source_tree)) {
        "MultiScholaR:::runProtDiaQcWorkerCommand"
    } else {
        "worker_env$runProtDiaQcWorkerCommand"
    }
    paste(
        "args <- commandArgs(trailingOnly = TRUE)",
        load,
        paste0(call, "(args[[2L]], args[[3L]])"),
        sep = "; "
    )
}

runProtDiaQcProcess <- function(spec, timeout_ms = 900000L) {
    spec <- validateProtDiaQcWorkerSpec(spec)
    spec_path <- tempfile("prot-dia-qc-spec-", fileext = ".rds")
    result_path <- tempfile("prot-dia-qc-result-", fileext = ".rds")
    on.exit(unlink(c(spec_path, result_path), force = TRUE), add = TRUE)
    saveRDS(spec, spec_path, version = 3L)
    namespace_path <- getNamespaceInfo(asNamespace("MultiScholaR"), "path")
    source_tree <- file.exists(file.path(
        namespace_path,
        "R",
        "mod_prot_dia_qc_worker_helpers.R"
    ))
    process <- processx::run(
        file.path(R.home("bin"), "Rscript"),
        c(
            "--vanilla", "-e", protDiaQcWorkerExpression(source_tree),
            namespace_path, spec_path, result_path
        ),
        error_on_status = FALSE,
        timeout = as.numeric(timeout_ms) / 1000,
        echo = FALSE
    )
    result <- if (file.exists(result_path)) readRDS(result_path) else NULL
    valid <- is.list(result) && isTRUE(result$ok) &&
        identical(as.integer(process$status), 0L) &&
        !isTRUE(result$complete_payload_returned)
    if (!isTRUE(valid)) {
        protDiaArtifactAbort(
            paste(
                "DIA-NN QC worker failed:",
                result$error_message %||% "missing worker result"
            ),
            "multischolar_prot_dia_qc_worker_failed"
        )
    }
    artifactResourceDataOnly(result, "DIA-NN QC worker result")
    result
}

applyProtDiaQcWorkerResult <- function(workflow_data, worker) {
    previous <- workflow_data$state_manager
    workflow_data$state_manager <- newProtDiaSettledStateManager(
        workflow_data,
        state_export = worker$state_export
    )
    if (inherits(previous, "ArtifactWorkflowState")) {
        try(previous$close(), silent = TRUE)
    }
    workflow_data$qc_params <- worker$qc_params
    recordProtDiaArtifactResult(workflow_data, worker$action_id, list(
        enabled = TRUE,
        ok = TRUE,
        committed = TRUE,
        stage_id = worker$action_id,
        state_name = worker$current_state,
        worker_pid = worker$worker_pid,
        complete_payload_returned = FALSE
    ))
    invisible(worker$result)
}

protDiaQcPlotGrid <- function(raw_png) {
    path <- tempfile("prot-dia-qc-parent-", fileext = ".png")
    on.exit(unlink(path, force = TRUE), add = TRUE)
    writeBin(raw_png, path)
    grid::rasterGrob(png::readPNG(path), interpolate = TRUE)
}
