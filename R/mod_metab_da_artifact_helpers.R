# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

.METAB_DA_RUN_SCHEMA <- "multischolar.metabolomics_da_run"
.METAB_DA_RUN_VERSION <- 1L
.METAB_DA_CURRENT_SCHEMA <- "multischolar.metabolomics_da_current"
.METAB_DA_CURRENT_VERSION <- 1L

metabDaArtifactMode <- function(kind = c("persistence", "queries")) {
    kind <- match.arg(kind)
    option <- paste0("multischolar.metabolomics.da_", kind)
    match.arg(getOption(option, "enabled"), c("enabled", "disabled"))
}

metabDaArtifactEligible <- function(workflow_data, kind = "persistence") {
    manager <- workflow_data$state_manager
    context <- workflow_data$workflow_context
    identical(metabDaArtifactMode(kind), "enabled") &&
        inherits(manager, "ArtifactWorkflowState") &&
        inherits(context, "WorkflowContext") && context$isBound() &&
        identical(context$getStorageDecision()$effective_backend, "artifact") &&
        identical(workflowStateType(manager), "metabolomics_standard")
}

metabDaPaths <- function(context, run_id = NULL) {
    root <- artifactNormalizeRelativePath(file.path(
        context$getPaths()$relative_paths$workflow_state_root,
        "differential_abundance"
    ))
    paths <- list(
        root = root,
        current = artifactNormalizeRelativePath(file.path(root, "current.json"))
    )
    if (!is.null(run_id)) {
        artifactValidatePathComponent(run_id, "metab_da_run_id")
        paths$run_root <- artifactNormalizeRelativePath(file.path(
            root,
            "runs",
            run_id
        ))
        paths$run_manifest <- artifactNormalizeRelativePath(file.path(
            paths$run_root,
            "manifest.json"
        ))
    }
    paths
}

metabDaStore <- function(context) {
    newArtifactStore(
        context$getPaths(),
        context$getIdentity()$project_id
    )
}

metabDaBaseColumns <- function() {
    c(
        "metabolite_id", "metabolite_name", "assay", "comparison",
        "friendly_name", "logFC", "raw_pvalue", "fdr_qvalue",
        "significant", "numerator", "denominator"
    )
}

metabDaValidateTable <- function(table, object, contrasts) {
    required <- metabDaBaseColumns()
    valid <- is.data.frame(table) && all(required %in% names(table)) &&
        !anyNA(table$assay) && !anyNA(table$comparison) &&
        !anyNA(table$metabolite_id) &&
        identical(sort(unique(table$assay)), sort(names(object@metabolite_data))) &&
        identical(sort(unique(table$comparison)), sort(contrasts))
    if (!isTRUE(valid)) {
        workflowSessionAbort(
            "metabolomics DA table is incomplete across assays or contrasts",
            "multischolar_incomplete_metab_da_run"
        )
    }
    keys <- artifactWorkflowStateEntityKeys(
        table,
        c("assay", "comparison", "metabolite_id")
    )
    if (length(keys) != nrow(table)) {
        workflowSessionAbort(
            "metabolomics DA query keys are ambiguous",
            "multischolar_ambiguous_metab_da_keys"
        )
    }
    invisible(table)
}

metabDaQuerySpecification <- function(table) {
    projections <- names(table)
    newArtifactQueryPageSpecification(
        query_id = "metabolomics.da.results.v1",
        state_role = "metab_da_results",
        projections = projections,
        filters = list(
            assay = list(
                column = "assay",
                type = "character",
                operators = c("equal", "in")
            ),
            comparison = list(
                column = "comparison",
                type = "character",
                operators = c("equal", "in")
            ),
            metabolite_id = list(
                column = "metabolite_id",
                type = "character",
                operators = c("equal", "in", "contains")
            ),
            metabolite_name = list(
                column = "metabolite_name",
                type = "character",
                operators = c("equal", "contains")
            ),
            fdr_qvalue = list(
                column = "fdr_qvalue",
                type = "double",
                operators = c("lt", "lte", "between")
            ),
            logFC = list(
                column = "logFC",
                type = "double",
                operators = c("gt", "gte", "lt", "lte", "between")
            ),
            significant = list(
                column = "significant",
                type = "character",
                operators = c("equal", "in")
            )
        ),
        sorts = list(
            fdr = list(
                column = "fdr_qvalue",
                transform = "identity",
                directions = c("asc", "desc")
            ),
            effect = list(
                column = "logFC",
                transform = "absolute",
                directions = c("asc", "desc")
            ),
            feature = list(
                column = "metabolite_id",
                transform = "identity",
                directions = c("asc", "desc")
            )
        ),
        default_sort = list(sort_id = "fdr", direction = "asc"),
        max_rows = 5000L,
        max_bytes = 64L * 1024L * 1024L
    )
}

metabDaRunId <- function(object, contrasts, parameters, results) {
    paste0("metab-da-", substr(artifactSemanticDigest(list(
        state_digest = metabNormSafeDigest(object),
        contrasts = contrasts,
        parameters = parameters,
        results = results$da_metabolites_long
    )), 1L, 24L))
}

metabDaManifestDigest <- function(value) {
    candidate <- value
    candidate$manifest_digest <- NULL
    workflowSessionContentDigest(candidate)
}

validateMetabDaManifest <- function(manifest) {
    required <- c(
        "schema", "schema_version", "run_id", "project_id", "workflow_id",
        "state_generation_id", "state_semantic_digest", "assays",
        "contrasts", "parameters", "software", "model_digest",
        "model_ref", "design_digest", "annotation_digest", "result_ref",
        "query_specification", "created_at", "manifest_digest"
    )
    valid <- is.list(manifest) && identical(names(manifest), required) &&
        identical(manifest$schema, .METAB_DA_RUN_SCHEMA) &&
        identical(
            workflowStateVersionValue(manifest$schema_version),
            .METAB_DA_RUN_VERSION
        ) && workflowCapabilityScalarString(manifest$run_id) &&
        workflowCapabilityScalarString(manifest$project_id) &&
        workflowCapabilityScalarString(manifest$workflow_id) &&
        workflowCapabilityScalarString(manifest$state_generation_id) &&
        is.character(unlist(manifest$assays, use.names = FALSE)) &&
        is.character(unlist(manifest$contrasts, use.names = FALSE)) &&
        is.list(manifest$parameters) && is.list(manifest$software) &&
        is.list(manifest$model_ref) &&
        is.list(manifest$result_ref) &&
        is.list(manifest$query_specification) &&
        artifactRefValidUtc(manifest$created_at)
    if (!isTRUE(valid)) {
        workflowSessionAbort(
            "metabolomics DA manifest is malformed",
            "multischolar_invalid_metab_da_manifest"
        )
    }
    for (field in c(
        "state_semantic_digest", "model_digest", "design_digest",
        "annotation_digest", "manifest_digest"
    )) {
        artifactRefValidateDigest(manifest[[field]], field)
    }
    artifactStoreNormalizeRef(manifest$result_ref)
    if (!identical(
        names(manifest$model_ref),
        c("relative_path", "byte_digest", "semantic_digest")
    )) {
        workflowSessionAbort(
            "metabolomics DA model ref is malformed",
            "multischolar_invalid_metab_da_manifest"
        )
    }
    artifactNormalizeRelativePath(manifest$model_ref$relative_path)
    artifactRefValidateDigest(manifest$model_ref$byte_digest, "model_byte_digest")
    artifactRefValidateDigest(
        manifest$model_ref$semantic_digest,
        "model_semantic_digest"
    )
    validateArtifactQueryPageSpecification(manifest$query_specification)
    if (!identical(manifest$manifest_digest, metabDaManifestDigest(manifest))) {
        workflowSessionAbort(
            "metabolomics DA manifest digest differs",
            "multischolar_metab_da_manifest_digest_mismatch"
        )
    }
    manifest
}

metabDaWriteModel <- function(store, relative_path, model, failure_injector) {
    path <- artifactStoreResolveFile(store, relative_path)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    temporary <- paste0(path, ".", artifactOpaqueId("tmp"), ".tmp")
    on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
    saveRDS(model, temporary, version = 3)
    check <- readRDS(temporary)
    if (!identical(check, model)) {
        workflowSessionAbort(
            "metabolomics DA model round trip changed content",
            "multischolar_metab_da_model_mismatch"
        )
    }
    artifactStoreInvokeFailure(
        failure_injector,
        "before_metab_da_model_publication",
        list(relative_path = relative_path)
    )
    if (!isTRUE(file.rename(temporary, path))) {
        workflowSessionAbort(
            "metabolomics DA model publication failed",
            "multischolar_metab_da_publish_failed"
        )
    }
    list(
        relative_path = relative_path,
        byte_digest = artifactByteDigest(path),
        semantic_digest = metabNormSafeDigest(model)
    )
}

metabDaWriteManifest <- function(
    manifest,
    path,
    replace = FALSE,
    failure_injector = NULL,
    failure_stage = "before_metab_da_manifest_publication"
) {
    validateMetabDaManifest(manifest)
    if (!isTRUE(replace) && file.exists(path)) {
        workflowSessionAbort(
            "immutable metabolomics DA manifest exists",
            "multischolar_metab_da_manifest_exists"
        )
    }
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    temporary <- paste0(path, ".", artifactOpaqueId("tmp"), ".tmp")
    on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
    writeLines(workflowSessionJson(unclass(manifest)), temporary)
    validateMetabDaManifest(
        jsonlite::read_json(temporary, simplifyVector = FALSE)
    )
    artifactStoreInvokeFailure(
        failure_injector,
        failure_stage,
        list(path = path)
    )
    if (!isTRUE(file.rename(temporary, path))) {
        workflowSessionAbort(
            "metabolomics DA manifest publication failed",
            "multischolar_metab_da_publish_failed"
        )
    }
    invisible(path)
}

persistMetabDaArtifacts <- function(
    workflow_data,
    results,
    contrasts_tbl,
    parameters,
    failure_injector = NULL
) {
    if (!metabDaArtifactEligible(workflow_data, "persistence")) {
        return(list(enabled = FALSE, ok = TRUE, reason = "artifact_da_disabled"))
    }
    object <- results$theObject
    table <- results$da_metabolites_long
    contrasts <- if ("contrasts" %in% names(contrasts_tbl)) {
        as.character(contrasts_tbl$contrasts)
    } else {
        as.character(contrasts_tbl$contrast_string)
    }
    metabDaValidateTable(table, object, contrasts)
    context <- workflow_data$workflow_context
    manager <- workflow_data$state_manager
    store <- metabDaStore(context)
    state_manifest <- artifactWorkflowStateReadManifest(
        store,
        manager$states[[manager$getCurrentStateName()]]$manifest_relative_path
    )
    run_id <- metabDaRunId(object, contrasts, parameters, results)
    paths <- metabDaPaths(context, run_id)
    encoded <- encodeArtifactTable(
        table,
        stable_key = c("assay", "comparison", "metabolite_id"),
        owner = "metabolomics.da.results"
    )
    ref <- artifactStoreWriteParquet(
        store,
        encoded,
        logical_key = list(
            project_id = context$getIdentity()$project_id,
            omic_type = "metabolomics",
            workflow_slug = context$getIdentity()$workflow_slug,
            stage_id = "differential_abundance",
            state_role = "metab_da_results",
            generation_id = run_id
        ),
        failure_injector = failure_injector
    )
    specification <- metabDaQuerySpecification(table)
    model_relative_path <- artifactNormalizeRelativePath(file.path(
        paths$run_root,
        "limma_models.rds"
    ))
    model_ref <- metabDaWriteModel(
        store,
        model_relative_path,
        results$contrasts_results,
        failure_injector
    )
    manifest <- list(
        schema = .METAB_DA_RUN_SCHEMA,
        schema_version = .METAB_DA_RUN_VERSION,
        run_id = run_id,
        project_id = context$getIdentity()$project_id,
        workflow_id = context$getIdentity()$workflow_id,
        state_generation_id = manager$getCurrentGenerationId(),
        state_semantic_digest = state_manifest$data$semantic_digest,
        assays = as.list(names(object@metabolite_data)),
        contrasts = as.list(contrasts),
        parameters = parameters,
        software = list(
            multischolar = metabQcPackageVersion(),
            limma = protDiaNormPackageVersion("limma"),
            r = as.character(getRversion())
        ),
        model_digest = metabNormSafeDigest(results$contrasts_results),
        model_ref = model_ref,
        design_digest = artifactSemanticDigest(object@design_matrix),
        annotation_digest = artifactSemanticDigest(lapply(
            object@metabolite_data,
            function(assay) assay[intersect(
                c(object@metabolite_id_column, object@annotation_id_column),
                names(assay)
            )]
        )),
        result_ref = ref,
        query_specification = specification,
        created_at = artifactRefUtcNow(),
        manifest_digest = NULL
    )
    manifest$manifest_digest <- metabDaManifestDigest(manifest)
    run_path <- artifactStoreResolveFile(store, paths$run_manifest)
    metabDaWriteManifest(manifest, run_path, failure_injector = failure_injector)
    current <- list(
        schema = .METAB_DA_CURRENT_SCHEMA,
        schema_version = .METAB_DA_CURRENT_VERSION,
        run_manifest_relative_path = paths$run_manifest,
        run_manifest_digest = manifest$manifest_digest,
        run_id = run_id,
        updated_at = artifactRefUtcNow()
    )
    current$manifest_digest <- workflowSessionContentDigest(current)
    current_path <- artifactStoreResolveFile(store, metabDaPaths(context)$current)
    temporary <- paste0(current_path, ".", artifactOpaqueId("tmp"), ".tmp")
    on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
    dir.create(dirname(current_path), recursive = TRUE, showWarnings = FALSE)
    writeLines(workflowSessionJson(current), temporary)
    artifactStoreInvokeFailure(
        failure_injector,
        "before_metab_da_current_publication",
        list(run_id = run_id)
    )
    if (!isTRUE(file.rename(temporary, current_path))) {
        workflowSessionAbort(
            "metabolomics DA current pointer publication failed",
            "multischolar_metab_da_publish_failed"
        )
    }
    result <- list(
        enabled = TRUE,
        ok = TRUE,
        committed = TRUE,
        run_id = run_id,
        manifest = manifest,
        result_ref = ref,
        query_specification = specification
    )
    recordArtifactStageResult(workflow_data, "differential_abundance", result)
    result
}

persistMetabDaArtifactsSafely <- function(...) {
    tryCatch(
        persistMetabDaArtifacts(...),
        error = function(error) list(
            enabled = TRUE,
            ok = FALSE,
            error = error,
            error_message = conditionMessage(error)
        )
    )
}

queryMetabDaPage <- function(
    workflow_data,
    projections = NULL,
    filters = list(),
    sort_id = NULL,
    direction = NULL,
    limit = 100L,
    cursor = NULL,
    resource_policy = NULL
) {
    result <- workflow_data$artifact_stage_results$differential_abundance
    if (!metabDaArtifactEligible(workflow_data, "queries") ||
        !is.list(result) || !isTRUE(result$committed)) {
        workflowSessionAbort(
            "metabolomics DA artifact query is unavailable",
            "multischolar_metab_da_query_unavailable"
        )
    }
    queryArtifactPage(
        metabDaStore(workflow_data$workflow_context),
        result$result_ref,
        result$query_specification,
        projections = projections,
        filters = filters,
        sort_id = sort_id,
        direction = direction,
        limit = limit,
        cursor = cursor,
        resource_policy = resource_policy
    )
}
