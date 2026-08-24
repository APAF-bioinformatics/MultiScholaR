# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

.PROT_DIA_SUMMARY_FINAL_STATES <- c(
    "correlation_filtered",
    "ruv_corrected",
    "protein_replicate_filtered",
    "imputed"
)
.PROT_DIA_SUMMARY_EXPORT_SCHEMA <- "multischolar.prot_dia.summary_export"
.PROT_DIA_SUMMARY_REPORT_SCHEMA <- "multischolar.prot_dia.report_dependencies"
.PROT_DIA_SUMMARY_EVIDENCE_VERSION <- 1L
.PROT_DIA_SUMMARY_MAX_REPORT_FILES <- 5000L

protDiaSummaryAbort <- function(message, class, ...) {
    rlang::abort(
        message,
        class = c(class, "multischolar_prot_dia_summary_error"),
        ...
    )
}

protDiaSummaryMode <- function(kind = c("report_reads", "final_export")) {
    kind <- match.arg(kind)
    option <- switch(kind,
        report_reads = "multischolar.prot_dia.summary_artifact_reads",
        final_export = "multischolar.prot_dia.summary_final_export"
    )
    match.arg(getOption(option, "enabled"), c("enabled", "disabled"))
}

protDiaSummaryArtifactEligible <- function(
    workflow_data,
    kind = c("report_reads", "final_export")
) {
    kind <- match.arg(kind)
    if (identical(protDiaSummaryMode(kind), "disabled") ||
        is.null(workflow_data)) {
        return(FALSE)
    }
    context <- tryCatch(
        workflow_data$workflow_context,
        error = \(error) NULL
    )
    manager <- tryCatch(
        workflow_data$state_manager,
        error = \(error) NULL
    )
    if (!inherits(context, "WorkflowContext") || !context$isBound() ||
        !inherits(manager, "ArtifactWorkflowState")) {
        return(FALSE)
    }
    identity <- context$getIdentity()
    decision <- context$getStorageDecision()
    identical(decision$effective_backend, "artifact") &&
        identical(workflowStateType(manager), "DIA") &&
        protDiaArtifactIdentityMatches(identity)
}

protDiaSummaryProjectPaths <- function(project_dirs, omic_type, context) {
    if (!is.list(project_dirs) || !workflowCapabilityScalarString(omic_type)) {
        protDiaSummaryAbort(
            "DIA-NN summary requires explicit project directories",
            "multischolar_missing_prot_dia_summary_dependency",
            dependency = "project_dirs"
        )
    }
    identity <- context$getIdentity()
    candidates <- unique(c(
        paste0(omic_type, "_", identity$omic_label),
        omic_type
    ))
    matches <- candidates[candidates %in% names(project_dirs)]
    if (length(matches) == 0L) {
        protDiaSummaryAbort(
            sprintf(
                "DIA-NN summary dependency 'project_dirs[%s]' is missing",
                omic_type
            ),
            "multischolar_missing_prot_dia_summary_dependency",
            dependency = "project_dirs"
        )
    }
    paths <- project_dirs[[matches[[1L]]]]
    required <- c("base_dir", "source_dir", "results_summary_dir")
    missing <- required[!vapply(required, function(name) {
        workflowCapabilityScalarString(paths[[name]])
    }, logical(1))]
    if (length(missing) > 0L) {
        protDiaSummaryAbort(
            sprintf(
                "DIA-NN summary project paths are missing: %s",
                paste(missing, collapse = ", ")
            ),
            "multischolar_missing_prot_dia_summary_dependency",
            dependency = paste0("project_dirs.", missing)
        )
    }
    project_root <- context$getProjectRoot()
    normalized <- lapply(paths[required], function(path) {
        normalizePath(path, winslash = "/", mustWork = FALSE)
    })
    outside <- required[!vapply(normalized, function(path) {
        artifactPathIsContained(path, project_root)
    }, logical(1))]
    if (length(outside) > 0L) {
        protDiaSummaryAbort(
            sprintf(
                "DIA-NN summary project path escapes its project root: %s",
                paste(outside, collapse = ", ")
            ),
            "multischolar_prot_dia_summary_path_escape",
            dependency = paste0("project_dirs.", outside)
        )
    }
    paths$qc_figures_dir <- paths$qc_figures_dir %||%
        file.path(paths$results_summary_dir, "QC_figures")
    paths$publication_figures_dir <- paths$publication_figures_dir %||%
        file.path(paths$results_summary_dir, "Publication_figures")
    paths$publication_tables_dir <- paths$publication_tables_dir %||%
        file.path(paths$results_summary_dir, "Publication_tables")
    paths
}

protDiaSummaryFinalStateName <- function(state_manager) {
    history <- workflowStateHistory(state_manager)
    state_name <- purrr::detect(
        .PROT_DIA_SUMMARY_FINAL_STATES,
        \(candidate) candidate %in% history
    )
    if (is.null(state_name)) {
        protDiaSummaryAbort(
            sprintf(
                "DIA-NN summary dependency 'final_s4_state' is missing; expected one of: %s",
                paste(.PROT_DIA_SUMMARY_FINAL_STATES, collapse = ", ")
            ),
            "multischolar_missing_prot_dia_summary_dependency",
            dependency = "final_s4_state"
        )
    }
    state_name
}

protDiaSummaryRequiredValue <- function(workflow_data, name) {
    value <- tryCatch(workflow_data[[name]], error = \(error) NULL)
    if (is.null(value)) {
        protDiaSummaryAbort(
            sprintf("DIA-NN summary dependency '%s' is missing", name),
            "multischolar_missing_prot_dia_summary_dependency",
            dependency = name
        )
    }
    value
}

protDiaSummaryRefEvidence <- function(refs) {
    lapply(refs, function(ref) {
        ref <- artifactStoreNormalizeRef(ref)
        list(
            artifact_id = ref$artifact_id,
            generation_id = ref$logical_key$generation_id,
            state_role = ref$logical_key$state_role,
            semantic_digest = ref$hash_policy$semantic$digest,
            byte_digest = ref$hash_policy$byte$digest
        )
    })
}

#' Resolve the exact summary workflow descriptor contract
#'
#' @param workflow_data Mutable proteomics workflow state.
#'
#' @return Descriptor identifier, version, and digest.
#' @noRd
protDiaSummaryDescriptorContract <- function(workflow_data) {
    context <- workflow_data$workflow_context
    descriptor <- findArtifactWorkflowDescriptor(
        context$getIdentity(),
        artifactWorkflowDescriptorCatalogue()
    )
    artifactStageDescriptorContract(descriptor)
}

protDiaSummaryDaAudit <- function(workflow_data) {
    index <- tryCatch(
        workflow_data$da_analysis_results_list,
        error = \(error) NULL
    )
    if (!isProtDiaDaArtifactIndex(index)) return(NULL)
    contrasts <- lapply(index$contrasts, function(entry) {
        list(
            contrast_id = entry$contrast_id,
            contrast_name = entry$contrast_name,
            full_format = entry$full_format,
            friendly_name = entry$friendly_name,
            manifest_digest = entry$manifest_digest,
            summary = entry$summary
        )
    })
    list(
        run_id = index$run_id,
        source_generation_id = index$source_generation_id,
        manifest_digest = index$manifest_digest,
        contrasts = contrasts
    )
}

protDiaSummaryEnrichmentAudit <- function(workflow_data) {
    index <- tryCatch(
        workflow_data$enrichment_artifact_index,
        error = \(error) NULL
    )
    if (!isProtDiaEnrichArtifactIndex(index)) return(NULL)
    requests <- lapply(index$requests, function(request) {
        list(
            request_id = request$request_id,
            request_digest = request$request_digest,
            backend = request$backend,
            contrast = request$contrast,
            direction = request$direction,
            status = request$status,
            execution_state = request$execution_state %||%
                tryCatch(request$service$source, error = \(error) NULL),
            attempts = request$attempts %||% request$attempt,
            response = request$response,
            legacy_service = request$service,
            legacy_software = request$software
        )
    })
    products <- lapply(index$products, function(product) {
        list(
            name = product$name,
            byte_digest = product$byte_digest,
            bytes = product$bytes
        )
    })
    list(
        run_id = index$run_id,
        source = index$source,
        parameters = index$parameters,
        software = index$software,
        manifest_digest = index$manifest_digest,
        requests = requests,
        products = products
    )
}

protDiaSummaryAnalysisAudit <- function(workflow_data) {
    list(
        da = protDiaSummaryDaAudit(workflow_data),
        enrichment = protDiaSummaryEnrichmentAudit(workflow_data),
        da_ui_params = tryCatch(
            workflow_data$da_ui_params,
            error = \(error) NULL
        ),
        enrichment_ui_params = tryCatch(
            workflow_data$enrichment_ui_params,
            error = \(error) NULL
        ),
        tab_status = tryCatch(
            workflow_data$tab_status[c(
                "differential_expression",
                "differential_abundance",
                "enrichment_analysis"
            )],
            error = \(error) NULL
        )
    )
}

protDiaSummarySessionAudit <- function(
    workflow_data,
    eligibility_fn = protDiaSummaryArtifactEligible
) {
    if (!eligibility_fn(workflow_data, "report_reads")) {
        return(NULL)
    }
    context <- tryCatch(
        workflow_data$workflow_context,
        error = \(error) NULL
    )
    manager <- tryCatch(
        workflow_data$state_manager,
        error = \(error) NULL
    )
    state_name <- protDiaSummaryFinalStateName(manager)
    metadata <- manager$getStateMetadata(state_name)
    design_matrix <- protDiaSummaryRequiredValue(workflow_data, "design_matrix")
    contrasts_tbl <- protDiaSummaryRequiredValue(workflow_data, "contrasts_tbl")
    audit <- list(
        identity = context$getIdentity()[c(
            "project_id", "workflow_id", "omic_type", "omic_label",
            "workflow_slug"
        )],
        descriptor_contract = protDiaSummaryDescriptorContract(workflow_data),
        final_state = list(
            state_name = state_name,
            generation_id = metadata$generation_id,
            s4_class = metadata$s4_class,
            artifact_refs = protDiaSummaryRefEvidence(metadata$artifact_refs),
            audit_metadata = metadata$audit_metadata
        ),
        scientific = protDiaSummaryScientificSummary(
            NULL,
            design_matrix,
            contrasts_tbl
        ),
        analysis = protDiaSummaryAnalysisAudit(workflow_data)
    )
    artifactWorkflowStateAssertSafeMetadata(
        audit,
        "DIA-NN summary session audit"
    )
    audit
}

protDiaSummaryValidateFinalObject <- function(object, design_matrix) {
    valid <- methods::is(object, "ProteinQuantitativeData") &&
        identical(methods::validObject(object, test = TRUE), TRUE) &&
        identical(object@design_matrix, design_matrix) &&
        nrow(object@protein_quant_table) > 0L &&
        nrow(object@protein_id_table) > 0L
    if (!isTRUE(valid)) {
        protDiaSummaryAbort(
            "DIA-NN final S4 export differs from its declared scientific dependencies",
            "multischolar_inexact_prot_dia_summary_export",
            dependency = "final_s4_object"
        )
    }
    invisible(TRUE)
}

protDiaSummaryAuditRecordId <- function(object, name) {
    if (is.null(object) || !methods::is(object, "ProteinQuantitativeData")) {
        return(NULL)
    }
    audit <- object@args[[name]]
    record_id <- audit$current_record_id
    if (workflowCapabilityScalarString(record_id)) record_id else NULL
}

protDiaSummaryScientificSummary <- function(
    final_object,
    design_matrix,
    contrasts_tbl
) {
    summary <- list(
        design_rows = as.integer(nrow(design_matrix)),
        design_columns = as.integer(ncol(design_matrix)),
        design_digest = artifactSemanticDigest(design_matrix),
        contrasts_digest = artifactSemanticDigest(contrasts_tbl)
    )
    if (is.null(final_object)) return(summary)
    c(summary, list(
        s4_class = class(final_object)[[1L]],
        protein_quant_rows = as.integer(nrow(final_object@protein_quant_table)),
        protein_quant_columns = as.integer(ncol(final_object@protein_quant_table)),
        protein_quant_digest = artifactSemanticDigest(
            final_object@protein_quant_table
        ),
        protein_id_rows = as.integer(nrow(final_object@protein_id_table)),
        protein_id_columns = as.integer(ncol(final_object@protein_id_table)),
        protein_id_digest = artifactSemanticDigest(final_object@protein_id_table),
        object_design_digest = artifactSemanticDigest(final_object@design_matrix),
        args_groups = names(final_object@args),
        args_digest = artifactSemanticDigest(final_object@args),
        peptide_audit_record_id = protDiaSummaryAuditRecordId(
            final_object,
            "peptide_qc_audit"
        ),
        protein_audit_record_id = protDiaSummaryAuditRecordId(
            final_object,
            "protein_qc_audit"
        ),
        normalization_audit_record_id = protDiaSummaryAuditRecordId(
            final_object,
            "normalization_artifact_audit"
        )
    ))
}

prepareProtDiaSummaryDependencies <- function(
    workflow_data,
    project_dirs,
    omic_type = "proteomics",
    kind = c("report_reads", "final_export"),
    hydrate_final = identical(match.arg(kind), "final_export"),
    eligibility_fn = protDiaSummaryArtifactEligible
) {
    kind <- match.arg(kind)
    if (!eligibility_fn(workflow_data, kind)) return(NULL)
    context <- workflow_data$workflow_context
    manager <- workflow_data$state_manager
    paths <- protDiaSummaryProjectPaths(project_dirs, omic_type, context)
    state_name <- protDiaSummaryFinalStateName(manager)
    state_metadata <- manager$getStateMetadata(state_name)
    design_matrix <- protDiaSummaryRequiredValue(workflow_data, "design_matrix")
    contrasts_tbl <- protDiaSummaryRequiredValue(workflow_data, "contrasts_tbl")
    final_object <- NULL
    if (isTRUE(hydrate_final)) {
        final_object <- manager$getState(state_name)
        protDiaSummaryValidateFinalObject(final_object, design_matrix)
    }
    config_list <- state_metadata$config
    if (!is.list(config_list)) {
        protDiaSummaryAbort(
            "DIA-NN summary dependency 'config_list' is missing from the final generation",
            "multischolar_missing_prot_dia_summary_dependency",
            dependency = "config_list"
        )
    }
    manifest <- list(
        identity = context$getIdentity()[c(
            "project_id", "workflow_id", "omic_type", "omic_label",
            "workflow_slug"
        )],
        descriptor_contract = protDiaSummaryDescriptorContract(workflow_data),
        final_state = list(
            state_name = state_name,
            generation_id = state_metadata$generation_id,
            s4_class = state_metadata$s4_class,
            artifact_refs = protDiaSummaryRefEvidence(
                state_metadata$artifact_refs
            ),
            audit_metadata = state_metadata$audit_metadata
        ),
        scientific = protDiaSummaryScientificSummary(
            final_object,
            design_matrix,
            contrasts_tbl
        ),
        analysis = protDiaSummaryAnalysisAudit(workflow_data)
    )
    artifactWorkflowStateAssertSafeMetadata(
        manifest,
        "DIA-NN summary dependency manifest"
    )
    structure(
        list(
            enabled = TRUE,
            kind = kind,
            projectRoot = context$getProjectRoot(),
            projectPaths = paths,
            finalS4Object = final_object,
            designMatrix = design_matrix,
            contrastsTbl = contrasts_tbl,
            configList = config_list,
            organismName = tryCatch(
                workflow_data$organism_name,
                error = \(error) NULL
            ),
            taxonId = tryCatch(
                workflow_data$taxon_id,
                error = \(error) NULL
            ),
            stateManager = manager,
            manifest = manifest
        ),
        class = c("MultiScholaRProtDiaSummaryDependencies", "list")
    )
}

releaseProtDiaSummaryDependencies <- function(dependencies) {
    if (!inherits(dependencies, "MultiScholaRProtDiaSummaryDependencies")) {
        return(invisible(FALSE))
    }
    dependencies$finalS4Object <- NULL
    workflowStateReleaseHydrationCache(dependencies$stateManager)
}

protDiaSummaryPathEvidence <- function(
    path,
    project_root,
    required,
    dependency
) {
    available <- file.exists(path) || dir.exists(path)
    if (!available && isTRUE(required)) {
        protDiaSummaryAbort(
            sprintf(
                "DIA-NN report dependency '%s' is missing",
                dependency
            ),
            "multischolar_missing_prot_dia_summary_dependency",
            dependency = dependency
        )
    }
    relative_path <- workflowSessionProjectRelativePath(path, project_root)
    if (!available) {
        return(list(
            relative_path = relative_path,
            required = required,
            available = FALSE,
            kind = "missing",
            files = 0L,
            bytes = 0,
            digest = NULL
        ))
    }
    if (!dir.exists(path)) {
        return(list(
            relative_path = relative_path,
            required = required,
            available = TRUE,
            kind = "file",
            files = 1L,
            bytes = unname(as.numeric(file.info(path)$size)),
            digest = artifactByteDigest(path)
        ))
    }
    files <- list.files(
        path,
        recursive = TRUE,
        all.files = TRUE,
        no.. = TRUE,
        full.names = TRUE
    )
    files <- files[!dir.exists(files)]
    if (length(files) > .PROT_DIA_SUMMARY_MAX_REPORT_FILES) {
        protDiaSummaryAbort(
            sprintf(
                "DIA-NN report dependency '%s' exceeds %d files",
                dependency,
                .PROT_DIA_SUMMARY_MAX_REPORT_FILES
            ),
            "multischolar_prot_dia_summary_dependency_too_large",
            dependency = dependency
        )
    }
    relative_files <- substring(files, nchar(path) + 2L)
    digests <- if (length(files) == 0L) {
        character()
    } else {
        vapply(files, artifactByteDigest, character(1))
    }
    bytes <- if (length(files) == 0L) {
        0
    } else {
        sum(as.numeric(file.info(files)$size))
    }
    list(
        relative_path = relative_path,
        required = required,
        available = TRUE,
        kind = "directory",
        files = as.integer(length(files)),
        bytes = unname(bytes),
        digest = artifactSemanticDigest(list(
            relative_files = relative_files,
            byte_digests = unname(digests)
        ))
    )
}

prepareProtDiaReportDependencies <- function(dependencies, template_path) {
    if (!inherits(dependencies, "MultiScholaRProtDiaSummaryDependencies")) {
        return(NULL)
    }
    paths <- dependencies$projectPaths
    root <- dependencies$projectRoot
    declared <- list(
        study_parameters = list(
            path = file.path(paths$source_dir, "study_parameters.txt"),
            required = TRUE
        ),
        design_matrix = list(
            path = file.path(paths$source_dir, "design_matrix.tab"),
            required = TRUE
        ),
        contrasts = list(
            path = file.path(paths$source_dir, "contrasts_tbl.tab"),
            required = FALSE
        ),
        template = list(path = template_path, required = TRUE),
        publication_tables = list(
            path = paths$publication_tables_dir,
            required = FALSE
        ),
        publication_figures = list(
            path = paths$publication_figures_dir,
            required = FALSE
        ),
        qc_figures = list(path = paths$qc_figures_dir, required = FALSE),
        pathway_enrichment = list(
            path = paths$pathway_dir %||% file.path(root, "pathway_enrichment"),
            required = FALSE
        )
    )
    evidence <- lapply(names(declared), function(name) {
        entry <- declared[[name]]
        protDiaSummaryPathEvidence(
            entry$path,
            root,
            entry$required,
            name
        )
    })
    names(evidence) <- names(declared)
    report_manifest <- list(
        schema = .PROT_DIA_SUMMARY_REPORT_SCHEMA,
        schema_version = .PROT_DIA_SUMMARY_EVIDENCE_VERSION,
        identity = dependencies$manifest$identity,
        descriptor_contract = dependencies$manifest$descriptor_contract,
        final_state = dependencies$manifest$final_state,
        scientific = dependencies$manifest$scientific,
        analysis = dependencies$manifest$analysis,
        dependencies = evidence,
        manifest_digest = NULL
    )
    report_manifest$manifest_digest <- artifactSemanticDigest(
        report_manifest
    )
    artifactWorkflowStateAssertSafeMetadata(
        report_manifest,
        "DIA-NN report dependency manifest"
    )
    dependencies$manifest <- report_manifest
    dependencies
}

protDiaSummaryValidateProduct <- function(object, expected) {
    valid <- identical(object, expected) &&
        methods::is(object, "ProteinQuantitativeData") &&
        identical(methods::validObject(object, test = TRUE), TRUE)
    if (!isTRUE(valid)) {
        protDiaSummaryAbort(
            "DIA-NN final S4 compatibility product failed exact validation",
            "multischolar_inexact_prot_dia_summary_export"
        )
    }
    invisible(TRUE)
}

protDiaSummaryWriteMetadata <- function(value, path) {
    artifactWorkflowStateAssertSafeMetadata(value, "DIA-NN summary metadata")
    encoded <- workflowSessionEncodeMetadata(value, "DIA-NN summary metadata")
    temporary <- file.path(
        dirname(path),
        paste0(".", basename(path), ".", artifactOpaqueId("tmp"), ".tmp")
    )
    on.exit(if (file.exists(temporary)) unlink(temporary, force = FALSE), add = TRUE)
    writeLines(encoded, temporary, useBytes = TRUE)
    decoded <- workflowSessionDecodeMetadata(
        paste(readLines(temporary, warn = FALSE), collapse = "\n"),
        "DIA-NN summary metadata"
    )
    if (!identical(decoded, value)) {
        protDiaSummaryAbort(
            "DIA-NN summary metadata changed during serialization",
            "multischolar_inexact_prot_dia_summary_metadata"
        )
    }
    protDiaSummaryPublishReplacement(
        temporary,
        path,
        "DIA-NN summary metadata"
    )
    invisible(path)
}

protDiaSummaryPublishReplacement <- function(temporary, path, owner) {
    backup <- NULL
    published <- FALSE
    if (file.exists(path) || dir.exists(path)) {
        backup <- file.path(
            dirname(path),
            paste0(".", basename(path), ".", artifactOpaqueId("old"), ".bak")
        )
        if (!isTRUE(file.rename(path, backup))) {
            protDiaSummaryAbort(
                sprintf("%s could not preserve its previous product", owner),
                "multischolar_prot_dia_summary_publish_failed",
                path = path
            )
        }
    }
    on.exit({
        if (!published && !is.null(backup) && file.exists(backup)) {
            try(file.rename(backup, path), silent = TRUE)
        }
    }, add = TRUE)
    if (!isTRUE(file.rename(temporary, path))) {
        protDiaSummaryAbort(
            sprintf("%s could not be published", owner),
            "multischolar_prot_dia_summary_publish_failed",
            path = path
        )
    }
    published <- TRUE
    if (!is.null(backup) && file.exists(backup)) {
        unlink(backup, force = FALSE)
    }
    invisible(path)
}

writeProtDiaSummaryFinalExport <- function(
    object,
    path,
    dependencies,
    save_rds_fn = saveRDS,
    read_rds_fn = readRDS
) {
    if (!inherits(dependencies, "MultiScholaRProtDiaSummaryDependencies")) {
        save_rds_fn(object, path)
        return(list(enabled = FALSE, path = path))
    }
    if (!artifactPathIsContained(path, dependencies$projectRoot)) {
        protDiaSummaryAbort(
            "DIA-NN final S4 export path escapes its project root",
            "multischolar_prot_dia_summary_path_escape"
        )
    }
    temporary <- file.path(
        dirname(path),
        paste0(".", basename(path), ".", artifactOpaqueId("tmp"), ".tmp")
    )
    on.exit(if (file.exists(temporary)) unlink(temporary, force = FALSE), add = TRUE)
    save_rds_fn(object, temporary)
    protDiaSummaryValidateProduct(read_rds_fn(temporary), object)
    protDiaSummaryPublishReplacement(
        temporary,
        path,
        "DIA-NN final S4 compatibility product"
    )
    product <- list(
        relative_path = workflowSessionProjectRelativePath(
            path,
            dependencies$projectRoot
        ),
        byte_digest = artifactByteDigest(path),
        bytes = unname(as.numeric(file.info(path)$size))
    )
    manifest <- list(
        schema = .PROT_DIA_SUMMARY_EXPORT_SCHEMA,
        schema_version = .PROT_DIA_SUMMARY_EVIDENCE_VERSION,
        identity = dependencies$manifest$identity,
        descriptor_contract = dependencies$manifest$descriptor_contract,
        final_state = dependencies$manifest$final_state,
        scientific = dependencies$manifest$scientific,
        analysis = dependencies$manifest$analysis,
        product = product,
        created_at = artifactRefUtcNow(),
        manifest_digest = NULL
    )
    manifest$manifest_digest <- artifactSemanticDigest(manifest)
    metadata_path <- paste0(path, ".artifact.json")
    protDiaSummaryWriteMetadata(manifest, metadata_path)
    list(
        enabled = TRUE,
        path = path,
        metadataPath = metadata_path,
        manifest = manifest
    )
}

recordProtDiaSummaryReportProduct <- function(dependencies, path) {
    if (!inherits(dependencies, "MultiScholaRProtDiaSummaryDependencies")) {
        return(NULL)
    }
    if (!file.exists(path) || dir.exists(path) ||
        !artifactPathIsContained(path, dependencies$projectRoot)) {
        protDiaSummaryAbort(
            "DIA-NN rendered report product is missing or outside its project",
            "multischolar_prot_dia_summary_report_product_mismatch"
        )
    }
    manifest <- dependencies$manifest
    manifest$report_product <- list(
        relative_path = workflowSessionProjectRelativePath(
            path,
            dependencies$projectRoot
        ),
        byte_digest = artifactByteDigest(path),
        bytes = unname(as.numeric(file.info(path)$size))
    )
    manifest$manifest_digest <- NULL
    manifest$manifest_digest <- artifactSemanticDigest(manifest)
    metadata_path <- paste0(path, ".artifact.json")
    protDiaSummaryWriteMetadata(manifest, metadata_path)
    list(path = path, metadataPath = metadata_path, manifest = manifest)
}
