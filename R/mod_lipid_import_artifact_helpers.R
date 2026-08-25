# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

lipidArtifactAbort <- function(message, class, ...) {
    rlang::abort(
        message,
        class = c(class, "multischolar_lipidomics_artifact_error"),
        ...
    )
}

lipidArtifactImportSpec <- function(format) {
    if (!identical(format, "lipidsearch")) return(NULL)
    list(
        capability_id = "lipidomics.lipidsearch.lipid.standard.v1",
        workflow_type = "lipidomics_standard",
        parser_id = "MultiScholaR::runLipidImportProcessing",
        parser_version = "1.0.0",
        format_id = "lipidsearch.direct_file",
        extensions = c("csv", "tsv", "txt")
    )
}

prepareLipidArtifactContext <- function(
    workflow_data,
    format = workflow_data$data_format,
    data_type = workflow_data$data_type,
    descriptor_catalogue = artifactWorkflowDescriptorCatalogue()
) {
    spec <- lipidArtifactImportSpec(format)
    if (is.null(spec) || !identical(data_type, "lipid")) {
        return(list(enabled = FALSE, reason = "not_supported_lipidomics_tuple"))
    }
    context <- workflow_data$workflow_context
    if (inherits(context, "WorkflowContext") &&
        !identical(context$getStaticIdentity()$omic_type, "lipidomics")) {
        return(list(enabled = FALSE, reason = "not_lipidomics_context"))
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
        lipidArtifactAbort(
            "lipidomics import resolved to the wrong exact descriptor",
            "multischolar_invalid_lipidomics_artifact_context"
        )
    }
    prepared
}

lipidArtifactCoordinatorOwned <- function(workflow_data) {
    artifactStageCoordinatorOwned(
        workflow_data,
        artifactLipidomicsWorkflowDescriptor()
    )
}

lipidArtifactAssayRoles <- function(assay_names, prefix = "assay") {
    roles <- sprintf("%s_%04d", prefix, seq_along(assay_names))
    stats::setNames(roles, assay_names)
}

lipidArtifactSourcePaths <- function(workflow_payload) {
    assay_names <- names(workflow_payload$assayList)
    explicit <- workflow_payload$assaySourceFiles
    if (!is.null(explicit)) {
        if (!identical(names(explicit), assay_names)) {
            lipidArtifactAbort(
                "lipidomics assay source names do not match assay order",
                "multischolar_invalid_lipidomics_source"
            )
        }
        return(unlist(explicit, use.names = TRUE))
    }
    sources <- workflow_payload$sourceFiles
    if (is.null(sources)) sources <- list()
    sources <- unlist(sources[!vapply(sources, is.null, logical(1))])
    if (length(sources) == 1L && length(assay_names) > 1L) {
        sources <- rep(sources, length(assay_names))
    }
    if (length(sources) != length(assay_names)) {
        lipidArtifactAbort(
            "each lipidomics assay requires one reviewed source binding",
            "multischolar_invalid_lipidomics_source"
        )
    }
    stats::setNames(unname(sources), assay_names)
}

lipidArtifactImportWorkflowPayload <- function(
    assay_list,
    column_mapping,
    data_format,
    sample_names_sanitized,
    assay1_file,
    assay2_file,
    processing_log
) {
    assay_names <- names(assay_list)
    source_files <- c(list(assay1_file), if (length(assay_names) > 1L) {
        list(assay2_file)
    } else {
        list()
    })
    if (length(source_files) != length(assay_names)) {
        lipidArtifactAbort(
            "lipidomics import source count differs from assay count",
            "multischolar_invalid_lipidomics_source"
        )
    }
    list(
        assayList = assay_list,
        sampleCols = column_mapping$sample_columns,
        sampleNamesSanitized = isTRUE(sample_names_sanitized),
        dataFormat = data_format,
        columnMapping = column_mapping,
        processingLog = processing_log,
        assaySourceFiles = stats::setNames(source_files, assay_names),
        sourceFiles = stats::setNames(source_files, assay_names)
    )
}

lipidArtifactValidateSource <- function(path, spec, assay_name) {
    valid <- workflowCapabilityScalarString(path) && file.exists(path) &&
        !dir.exists(path)
    extension <- if (isTRUE(valid)) tolower(tools::file_ext(path)) else ""
    if (!isTRUE(valid) || !extension %in% spec$extensions) {
        lipidArtifactAbort(
            sprintf("assay '%s' has no reviewed LipidSearch source", assay_name),
            "multischolar_unverified_lipidomics_source_encoding",
            assay_name = assay_name,
            source_extension = extension
        )
    }
    path
}

lipidArtifactAssayIdentity <- function(
    workflow_id,
    capability_id,
    assay_name,
    assay_order,
    source_digest,
    table_digest,
    mapping_digest
) {
    paste0("lipassay_", artifactSemanticDigest(list(
        workflow_id = workflow_id,
        capability_id = capability_id,
        assay_name = assay_name,
        assay_order = as.integer(assay_order),
        source_digest = source_digest,
        table_digest = table_digest,
        mapping_digest = mapping_digest
    )))
}

lipidArtifactImportManifest <- function(
    workflow_payload,
    spec,
    source_paths,
    identity
) {
    assay_names <- names(workflow_payload$assayList)
    roles <- lipidArtifactAssayRoles(assay_names)
    mapping <- workflow_payload$columnMapping
    mapping_digest <- artifactSemanticDigest(mapping)
    rows <- lapply(seq_along(assay_names), \(index) {
        assay_name <- assay_names[[index]]
        source_path <- lipidArtifactValidateSource(
            source_paths[[assay_name]],
            spec,
            assay_name
        )
        assay <- workflow_payload$assayList[[assay_name]]
        sample_columns <- mapping$sample_columns
        source_digest <- artifactByteDigest(source_path)
        table_digest <- artifactSemanticDigest(assay)
        data.frame(
            assay_name = assay_name,
            assay_order = as.integer(index),
            assay_identity = lipidArtifactAssayIdentity(
                identity$workflow_id,
                spec$capability_id,
                assay_name,
                index,
                source_digest,
                table_digest,
                mapping_digest
            ),
            artifact_role = roles[[assay_name]],
            source_digest = source_digest,
            source_size_bytes = unname(as.numeric(file.info(source_path)$size)),
            parser_id = spec$parser_id,
            parser_version = spec$parser_version,
            format_id = spec$format_id,
            mapping_digest = mapping_digest,
            table_digest = table_digest,
            feature_count = as.integer(nrow(assay)),
            sample_count = as.integer(length(sample_columns)),
            lipid_id_column = mapping$lipid_id_col,
            lipid_class_column = mapping$annotation_col %||% NA_character_,
            sample_columns_json = artifactStageParameterJson(sample_columns),
            column_schema_json = artifactStageParameterJson(names(assay)),
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, rows)
}

lipidArtifactImportTables <- function(workflow_payload, manifest) {
    assay_names <- manifest$assay_name
    assay_tables <- workflow_payload$assayList[assay_names]
    names(assay_tables) <- manifest$artifact_role
    c(assay_tables, list(
        assay_manifest = manifest,
        column_mapping = artifactStageMetadataTable(
            workflow_payload$columnMapping
        ),
        source_manifest = manifest[c(
            "assay_name", "assay_order", "source_digest",
            "source_size_bytes", "parser_id", "parser_version", "format_id"
        )]
    ))
}

lipidArtifactAggregateSource <- function(manifest) {
    unique_sources <- !duplicated(manifest$source_digest)
    list(
        source_role = "lipidomics_assay_set",
        source_uri = NULL,
        source_digest = artifactSemanticDigest(list(
            assay_names = manifest$assay_name,
            source_digests = manifest$source_digest
        )),
        source_size_bytes = sum(manifest$source_size_bytes[unique_sources]),
        parser_id = unique(manifest$parser_id)[[1L]],
        parser_version = unique(manifest$parser_version)[[1L]],
        format_id = unique(manifest$format_id)[[1L]],
        data_level = "lipid"
    )
}

persistLipidImportArtifacts <- function(
    workflow_data,
    workflow_payload,
    failure_injector = NULL,
    log_warn = logger::log_warn
) {
    runArtifactStageSafely(
        workflow_data,
        "import",
        \() {
            prepared <- prepareLipidArtifactContext(workflow_data)
            if (!isTRUE(prepared$enabled)) {
                return(list(
                    enabled = FALSE,
                    ok = TRUE,
                    stage_id = "import",
                    reason = prepared$reason
                ))
            }
            spec <- lipidArtifactImportSpec(workflow_data$data_format)
            exact_memory <- identical(
                workflow_data$data_tbl,
                workflow_payload$assayList
            ) && identical(
                workflow_data$column_mapping,
                workflow_payload$columnMapping
            )
            if (is.null(spec) || !isTRUE(exact_memory)) {
                lipidArtifactAbort(
                    "lipidomics artifact import differs from memory handoff",
                    "multischolar_inexact_lipidomics_import"
                )
            }
            sources <- lipidArtifactSourcePaths(workflow_payload)
            manifest <- lipidArtifactImportManifest(
                workflow_payload,
                spec,
                sources,
                prepared$context$getIdentity()
            )
            stage <- writeArtifactStage(
                prepared$context,
                prepared$descriptor,
                stage_id = "import",
                tables = lipidArtifactImportTables(workflow_payload, manifest),
                parameters = list(
                    capability_id = spec$capability_id,
                    assay_order = manifest$assay_name,
                    assay_roles = manifest$artifact_role,
                    assay_identities = manifest$assay_identity,
                    column_mapping = workflow_payload$columnMapping,
                    column_mapping_serialized =
                        artifactWorkflowStateSerializeMetadata(
                            workflow_payload$columnMapping,
                            "lipidomics import column mapping"
                        ),
                    readthrough_contract_version = 1L,
                    input_format = "lipidsearch",
                    data_level = "lipid",
                    sample_names_sanitized = isTRUE(
                        workflow_payload$sampleNamesSanitized
                    ),
                    evidence_boundary = paste(
                        "generated inputs certify scale only;",
                        "reviewed fixtures certify parser and scientific parity"
                    )
                ),
                source = lipidArtifactAggregateSource(manifest),
                failure_injector = failure_injector,
                abort_fn = lipidArtifactAbort
            )
            assay_refs <- stage$refs[manifest$artifact_role]
            names(assay_refs) <- manifest$assay_name
            c(
                list(enabled = TRUE, ok = TRUE, assay_refs = assay_refs),
                stage,
                list(assay_manifest = manifest)
            )
        },
        "lipidomics",
        log_warn
    )
}
