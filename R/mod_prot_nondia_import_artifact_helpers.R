# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

#' Abort non-DIA proteomics artifact persistence
#'
#' @param message Human-readable failure message.
#' @param class Specific condition class.
#' @param ... Additional fields passed to [rlang::abort()].
#'
#' @return This function does not return; it signals a typed error.
#' @noRd
protNonDiaArtifactAbort <- function(message, class, ...) {
    rlang::abort(
        message,
        class = c(class, "multischolar_prot_nondia_artifact_error"),
        ...
    )
}

#' Define exact non-DIA proteomics import contracts
#'
#' @return A named list of parser, format, and parameter specifications.
#' @noRd
protNonDiaArtifactImportSpecs <- function() {
    list(
        maxquant = list(
            workflow_type = "LFQ",
            capability_id = "proteomics.maxquant.protein.lfq.v1",
            parser_id = "MultiScholaR::importMaxQuantData",
            parser_version = "1.0.0",
            format_id = "maxquant.protein_groups.txt",
            extension = "txt",
            parser_parameters = list(
                use_lfq = TRUE,
                filter_contaminants = TRUE
            )
        ),
        fragpipe = list(
            workflow_type = "LFQ",
            capability_id = "proteomics.fragpipe.protein.lfq.v1",
            parser_id = "MultiScholaR::importFragPipeData",
            parser_version = "1.0.0",
            format_id = "fragpipe.combined_protein.tsv",
            extension = "tsv",
            parser_parameters = list(use_maxlfq = TRUE)
        ),
        pd_tmt = list(
            workflow_type = "TMT",
            capability_id = "proteomics.pd_tmt.protein.tmt.v1",
            parser_id = "MultiScholaR::importProteomeDiscovererTMTData",
            parser_version = "1.0.0",
            format_id = "proteome_discoverer.tmt.tsv",
            extension = "tsv",
            parser_parameters = list()
        )
    )
}

#' Resolve a non-DIA proteomics import contract
#'
#' @param format Candidate import format identifier.
#'
#' @return The matching import specification, or `NULL` when unsupported.
#' @noRd
protNonDiaArtifactImportSpec <- function(format) {
    if (!workflowCapabilityScalarString(format)) return(NULL)
    protNonDiaArtifactImportSpecs()[[format]]
}

#' Prepare an exact non-DIA proteomics artifact context
#'
#' @param workflow_data Mutable proteomics workflow state.
#' @param format Candidate import format.
#' @param data_type Candidate proteomics data level.
#' @param descriptor_catalogue Workflow descriptor catalogue.
#'
#' @return A prepared context result with an explicit enabled flag and reason.
#' @noRd
prepareProtNonDiaArtifactContext <- function(
    workflow_data,
    format = workflow_data$data_format,
    data_type = workflow_data$data_type,
    descriptor_catalogue = artifactWorkflowDescriptorCatalogue()
) {
    spec <- protNonDiaArtifactImportSpec(format)
    if (is.null(spec) || !identical(data_type, "protein")) {
        return(list(enabled = FALSE, reason = "not_supported_nondia_protein"))
    }
    context <- workflow_data$workflow_context
    if (inherits(context, "WorkflowContext") &&
        !identical(context$getStaticIdentity()$omic_type, "proteomics")) {
        return(list(enabled = FALSE, reason = "not_proteomics_context"))
    }
    prepared <- prepareArtifactStageContext(
        workflow_data,
        workflow_type = spec$workflow_type,
        input_format = format,
        data_level = data_type,
        descriptor_catalogue = descriptor_catalogue
    )
    if (isTRUE(prepared$enabled) &&
        !identical(prepared$descriptor$descriptor_id, spec$capability_id)) {
        protNonDiaArtifactAbort(
            "non-DIA import resolved to the wrong exact workflow descriptor",
            "multischolar_invalid_prot_nondia_artifact_context"
        )
    }
    prepared
}

#' Safely prepare a non-DIA proteomics artifact context
#'
#' @param workflow_data Mutable proteomics workflow state.
#' @param format Candidate import format.
#' @param data_type Candidate proteomics data level.
#' @param log_warn Warning logger used for additive artifact failures.
#'
#' @return The recorded context preparation result, invisibly.
#' @noRd
prepareProtNonDiaArtifactContextSafely <- function(
    workflow_data,
    format = workflow_data$data_format,
    data_type = workflow_data$data_type,
    log_warn = logger::log_warn
) {
    runArtifactStageSafely(
        workflow_data,
        "context",
        \() {
            result <- prepareProtNonDiaArtifactContext(
                workflow_data,
                format,
                data_type
            )
            result$ok <- TRUE
            result$stage_id <- "context"
            result
        },
        "non-DIA proteomics",
        log_warn
    )
}

#' Dispatch proteomics artifact context preparation by exact tuple
#'
#' @param workflow_data Mutable proteomics workflow state.
#' @param format Candidate import format.
#' @param data_type Candidate proteomics data level.
#' @param log_warn Warning logger used for additive artifact failures.
#'
#' @return The recorded DIA or non-DIA context preparation result, invisibly.
#' @noRd
prepareProtArtifactContextSafely <- function(
    workflow_data,
    format = workflow_data$data_format,
    data_type = workflow_data$data_type,
    log_warn = logger::log_warn
) {
    if (identical(format, "diann") && identical(data_type, "peptide")) {
        return(prepareProtDiaArtifactContextSafely(
            workflow_data,
            format,
            data_type,
            log_warn
        ))
    }
    prepareProtNonDiaArtifactContextSafely(
        workflow_data,
        format,
        data_type,
        log_warn
    )
}

#' Capture redacted provenance for a non-DIA import source
#'
#' @param source_path Existing reviewed import file.
#' @param format Exact supported import format.
#' @param retain_source_uri Whether to retain the source basename.
#'
#' @return A source registry record containing parser and byte-level provenance.
#' @noRd
protNonDiaArtifactSourceMetadata <- function(
    source_path,
    format,
    retain_source_uri = FALSE
) {
    spec <- protNonDiaArtifactImportSpec(format)
    valid_source <- artifactResourceScalarString(source_path) &&
        file.exists(source_path) && !dir.exists(source_path)
    if (is.null(spec) || !isTRUE(valid_source)) {
        protNonDiaArtifactAbort(
            "non-DIA proteomics artifact source must be one supported file",
            "multischolar_invalid_prot_nondia_source"
        )
    }
    extension <- tolower(tools::file_ext(source_path))
    if (!identical(extension, spec$extension)) {
        protNonDiaArtifactAbort(
            sprintf(
                "non-DIA tuple '%s' has no artifact proof for '.%s' input",
                spec$capability_id,
                extension
            ),
            "multischolar_unverified_prot_nondia_source_encoding",
            capability_id = spec$capability_id,
            source_extension = extension,
            remediation = sprintf(
                "Use the reviewed '.%s' encoding or remain memory-backed.",
                spec$extension
            )
        )
    }
    list(
        source_role = "search_results",
        source_uri = if (isTRUE(retain_source_uri)) basename(source_path) else NULL,
        source_digest = artifactByteDigest(source_path),
        source_size_bytes = unname(as.numeric(file.info(source_path)$size)),
        parser_id = spec$parser_id,
        parser_version = spec$parser_version,
        format_id = spec$format_id,
        data_level = "protein"
    )
}

#' Dual-write a completed non-DIA proteomics import
#'
#' @param workflow_data Mutable proteomics workflow state after memory import.
#' @param data_import_result Result returned by the reviewed scientific parser.
#' @param source_path Existing import source file.
#' @param sanitize_names Whether run names were sanitized by the memory workflow.
#' @param retain_source_uri Whether to retain the source basename.
#' @param failure_injector Optional artifact failure injector used by tests.
#' @param log_warn Warning logger used for additive artifact failures.
#'
#' @return The recorded import persistence result, invisibly.
#' @noRd
persistProtNonDiaImportArtifacts <- function(
    workflow_data,
    data_import_result,
    source_path,
    sanitize_names = FALSE,
    retain_source_uri = FALSE,
    failure_injector = NULL,
    log_warn = logger::log_warn
) {
    runArtifactStageSafely(
        workflow_data,
        "import",
        \() {
            format <- workflow_data$data_format
            spec <- protNonDiaArtifactImportSpec(format)
            prepared <- prepareProtNonDiaArtifactContext(workflow_data)
            if (!isTRUE(prepared$enabled)) {
                return(list(
                    enabled = FALSE,
                    ok = TRUE,
                    stage_id = "import",
                    reason = prepared$reason
                ))
            }
            if (is.null(spec) ||
                !identical(data_import_result$data_type, "protein") ||
                !identical(data_import_result$column_mapping,
                    workflow_data$column_mapping)) {
                protNonDiaArtifactAbort(
                    "non-DIA artifact import metadata differs from memory import",
                    "multischolar_inexact_prot_nondia_import_metadata"
                )
            }
            stage <- writeArtifactStage(
                prepared$context,
                prepared$descriptor,
                stage_id = "import",
                tables = list(canonical_data = workflow_data$data_tbl),
                parameters = list(
                    capability_id = spec$capability_id,
                    column_mapping = data_import_result$column_mapping,
                    column_mapping_serialized =
                        artifactWorkflowStateSerializeMetadata(
                            data_import_result$column_mapping,
                            "non-DIA import column mapping"
                        ),
                    readthrough_contract_version = 1L,
                    workflow_type = spec$workflow_type,
                    input_format = format,
                    data_level = "protein",
                    parser_parameters = spec$parser_parameters,
                    sanitize_names = isTRUE(sanitize_names)
                ),
                source = protNonDiaArtifactSourceMetadata(
                    source_path,
                    format,
                    retain_source_uri
                ),
                failure_injector = failure_injector,
                abort_fn = protNonDiaArtifactAbort
            )
            c(list(enabled = TRUE, ok = TRUE), stage)
        },
        "non-DIA proteomics",
        log_warn
    )
}

#' Dispatch proteomics import dual-write by exact workflow tuple
#'
#' @param workflow_data Mutable proteomics workflow state after memory import.
#' @param data_import_result Result returned by the reviewed scientific parser.
#' @param source_path Existing import source file.
#' @param use_precursor_norm Whether DIA-NN used normalized precursor abundance.
#' @param sanitize_names Whether run names were sanitized by the memory workflow.
#' @param retain_source_uri Whether to retain the source basename.
#' @param failure_injector Optional artifact failure injector used by tests.
#' @param pending_stage Optional independently verified DIA worker stage.
#' @param worker_attempted Whether the DIA worker path was attempted.
#' @param log_warn Warning logger used for additive artifact failures.
#'
#' @return The recorded DIA or non-DIA import persistence result, invisibly.
#' @noRd
persistProtImportArtifacts <- function(
    workflow_data,
    data_import_result,
    source_path,
    use_precursor_norm = TRUE,
    sanitize_names = FALSE,
    retain_source_uri = FALSE,
    failure_injector = NULL,
    pending_stage = NULL,
    worker_attempted = !is.null(pending_stage),
    log_warn = logger::log_warn
) {
    if (identical(workflow_data$data_format, "diann") &&
        identical(workflow_data$data_type, "peptide")) {
        return(persistProtDiaImportArtifacts(
            workflow_data,
            data_import_result,
            source_path,
            use_precursor_norm,
            sanitize_names,
            retain_source_uri,
            failure_injector,
            pending_stage,
            worker_attempted,
            log_warn
        ))
    }
    persistProtNonDiaImportArtifacts(
        workflow_data,
        data_import_result,
        source_path,
        sanitize_names,
        retain_source_uri,
        failure_injector,
        log_warn
    )
}
