# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

.PROT_DIA_DA_RUN_SCHEMA <- "multischolar.prot_dia_da_run"
.PROT_DIA_DA_RUN_VERSION <- 1L
.PROT_DIA_DA_CONTRAST_SCHEMA <- "multischolar.prot_dia_da_contrast"
.PROT_DIA_DA_CONTRAST_VERSION <- 1L
.PROT_DIA_DA_MODEL_SCHEMA <- "multischolar.prot_dia_da_model_rds"
.PROT_DIA_DA_MODEL_VERSION <- 1L
.PROT_DIA_DA_INDEX_SCHEMA <- "multischolar.prot_dia_da_index"
.PROT_DIA_DA_INDEX_VERSION <- 1L
.PROT_DIA_DA_MAX_QUERY_ROWS <- 250000L
.PROT_DIA_DA_MAX_QUERY_BYTES <- 256L * 1024L * 1024L

protDiaDaArtifactAbort <- function(message, class, ...) {
    rlang::abort(
        message,
        class = c(class, "multischolar_prot_dia_da_artifact_error"),
        ...
    )
}

protDiaDaArtifactMode <- function(kind = c("persistence", "queries")) {
    kind <- match.arg(kind)
    option <- paste0("multischolar.prot_dia.da_", kind)
    match.arg(getOption(option, "enabled"), c("enabled", "disabled"))
}

protDiaDaArtifactWorkflow <- function(workflow_data) {
    if (is.null(workflow_data)) return(FALSE)
    manager <- tryCatch(
        workflow_data$state_manager,
        error = function(error) NULL
    )
    context <- tryCatch(
        workflow_data$workflow_context,
        error = function(error) NULL
    )
    if (!inherits(manager, "ArtifactWorkflowState") ||
        !inherits(context, "WorkflowContext") || !context$isBound()) {
        return(FALSE)
    }
    identity <- context$getIdentity()
    decision <- context$getStorageDecision()
    identical(decision$effective_backend, "artifact") &&
        identical(workflowStateType(manager), "DIA") &&
        protDiaArtifactIdentityMatches(identity)
}

protDiaDaArtifactEligible <- function(workflow_data, kind = "persistence") {
    identical(protDiaDaArtifactMode(kind), "enabled") &&
        protDiaDaArtifactWorkflow(workflow_data)
}

protDiaDaSoftware <- function() {
    list(
        multischolar = protDiaNormPackageVersion("MultiScholaR"),
        limma = protDiaNormPackageVersion("limma"),
        qvalue = protDiaNormPackageVersion("qvalue"),
        r = as.character(getRversion())
    )
}

protDiaDaParameters <- function(
    formula_string,
    da_q_val_thresh,
    treat_lfc_cutoff,
    eBayes_trend = TRUE,
    eBayes_robust = TRUE
) {
    parameters <- list(
        formula_string = formula_string,
        da_q_val_thresh = as.double(da_q_val_thresh),
        treat_lfc_cutoff = as.double(treat_lfc_cutoff),
        eBayes_trend = isTRUE(eBayes_trend),
        eBayes_robust = isTRUE(eBayes_robust),
        qvalue_column = "fdr_qvalue",
        raw_pvalue_column = "raw_pvalue",
        fdr_fallback = "BH"
    )
    artifactWorkflowStateAssertSafeMetadata(parameters, "DIA-NN DA parameters")
    parameters
}

protDiaDaPaths <- function(context, run_id = NULL, contrast_id = NULL) {
    relative_root <- artifactNormalizeRelativePath(file.path(
        context$getPaths()$relative_paths$workflow_state_root,
        "differential_abundance"
    ))
    paths <- list(
        root = relative_root,
        current = artifactNormalizeRelativePath(file.path(
            relative_root, "current.json"
        ))
    )
    if (!is.null(run_id)) {
        artifactValidatePathComponent(run_id, "da_run_id")
        paths$run_root <- artifactNormalizeRelativePath(file.path(
            relative_root, "runs", run_id
        ))
        paths$run_manifest <- artifactNormalizeRelativePath(file.path(
            paths$run_root, "manifest.json"
        ))
    }
    if (!is.null(contrast_id)) {
        artifactValidatePathComponent(contrast_id, "da_contrast_id")
        paths$contrast_root <- artifactNormalizeRelativePath(file.path(
            paths$run_root, "contrasts", contrast_id
        ))
        paths$contrast_manifest <- artifactNormalizeRelativePath(file.path(
            paths$contrast_root, "manifest.json"
        ))
        paths$model <- artifactNormalizeRelativePath(file.path(
            paths$contrast_root, "limma_model.rds"
        ))
    }
    paths
}

protDiaDaStore <- function(context) {
    identity <- context$getIdentity()
    newArtifactStore(context$getPaths(), identity$project_id)
}

protDiaDaJsonDigest <- function(value) {
    workflowSessionContentDigest(value)
}

protDiaDaWriteJson <- function(
    value,
    path,
    validator,
    replace = FALSE,
    failure_injector = NULL,
    failure_stage = "before_da_json_publication"
) {
    validator(value)
    if (!isTRUE(replace) && (file.exists(path) || dir.exists(path))) {
        protDiaDaArtifactAbort(
            "immutable DIA-NN DA manifest already exists",
            "multischolar_prot_dia_da_manifest_exists",
            path = path
        )
    }
    parent <- dirname(path)
    if (!dir.exists(parent) &&
        !dir.create(parent, recursive = TRUE, showWarnings = FALSE)) {
        protDiaDaArtifactAbort(
            "DIA-NN DA manifest directory could not be created",
            "multischolar_prot_dia_da_write_failed"
        )
    }
    temporary <- file.path(
        parent,
        paste0(".", basename(path), ".", artifactOpaqueId("tmp"), ".tmp")
    )
    on.exit(if (file.exists(temporary)) unlink(temporary, force = FALSE), add = TRUE)
    connection <- file(temporary, open = "wb")
    tryCatch(
        writeBin(
            charToRaw(paste0(workflowSessionJson(unclass(value)), "\n")),
            connection
        ),
        finally = close(connection)
    )
    decoded <- jsonlite::read_json(temporary, simplifyVector = FALSE)
    validator(decoded)
    artifactStoreInvokeFailure(
        failure_injector,
        failure_stage,
        list(path = path)
    )
    if (!isTRUE(file.rename(temporary, path))) {
        protDiaDaArtifactAbort(
            "DIA-NN DA manifest could not be atomically published",
            "multischolar_prot_dia_da_publish_failed",
            path = path
        )
    }
    invisible(path)
}

protDiaDaReadJson <- function(path, validator) {
    if (!workflowCapabilityScalarString(path) || !file.exists(path) ||
        dir.exists(path)) {
        protDiaDaArtifactAbort(
            "DIA-NN DA manifest does not exist",
            "multischolar_missing_prot_dia_da_manifest",
            path = path
        )
    }
    value <- tryCatch(
        jsonlite::read_json(path, simplifyVector = FALSE),
        error = function(error) protDiaDaArtifactAbort(
            "DIA-NN DA manifest could not be parsed",
            "multischolar_malformed_prot_dia_da_manifest",
            parent = error
        )
    )
    validator(value)
}

protDiaDaStageEvidence <- function(workflow_data) {
    session <- workflow_data$artifact_session_manifest
    if (is.list(session) && is.list(session$stage_refs)) {
        return(session$stage_refs)
    }
    protDaSessionStageEvidence(workflow_data)
}

protDiaDaSourceBinding <- function(workflow_data) {
    manager <- workflow_data$state_manager
    context <- workflow_data$workflow_context
    snapshot <- workflowStateManifest(manager)
    current <- tail(snapshot$active_lineage, 1L)[[1L]]
    store <- protDiaDaStore(context)
    manifest <- artifactWorkflowStateReadManifest(
        store, current$manifest_relative_path
    )
    if (!identical(current$generation_id, snapshot$current_generation_id) ||
        !identical(current$logical_name, snapshot$current_state)) {
        protDiaDaArtifactAbort(
            "DIA-NN DA source lineage does not end at the current state",
            "multischolar_prot_dia_da_source_mismatch"
        )
    }
    list(
        generation_id = current$generation_id,
        state_name = current$logical_name,
        manifest_relative_path = current$manifest_relative_path,
        manifest_digest = manifest$manifest_digest,
        state_semantic_digest = manifest$data$semantic_digest,
        stage_refs = protDiaDaStageEvidence(workflow_data)
    )
}

protDiaDaRunId <- function(source, contrasts, parameters, software) {
    paste0("darun_", artifactSemanticDigest(list(
        source_generation_id = source$generation_id,
        source_state_digest = source$state_semantic_digest,
        stage_refs = source$stage_refs,
        contrasts = contrasts,
        parameters = parameters,
        software = software
    )))
}

protDiaDaRunContrasts <- function(contrasts) {
    data.frame(
        contrasts = as.character(contrasts$contrasts),
        full_format = as.character(contrasts$full_format),
        friendly_names = as.character(contrasts$friendly_names),
        stringsAsFactors = FALSE
    )
}

protDiaDaContrastMetadata <- function(contrasts, index) {
    list(
        contrast_name = as.character(contrasts$contrasts[[index]]),
        full_format = as.character(contrasts$full_format[[index]]),
        friendly_name = as.character(contrasts$friendly_names[[index]])
    )
}

protDiaDaContrastId <- function(metadata) {
    paste0("contrast_", substr(artifactSemanticDigest(metadata), 1L, 24L))
}

protDiaDaResultTables <- function(result) {
    tables <- result[vapply(result, is.data.frame, logical(1))]
    count_table <- result$num_sig_da_molecules_first_go$table
    if (is.data.frame(count_table)) {
        tables$num_sig_da_molecules <- count_table
    }
    tables
}

protDiaDaValidateLongTable <- function(table, protein_id_column, metadata) {
    required <- c(
        protein_id_column, "comparison", "log2FC", "raw_pvalue",
        "fdr_qvalue"
    )
    if (!is.data.frame(table) || !all(required %in% names(table))) {
        protDiaDaArtifactAbort(
            "DIA-NN DA long table is missing its scientific result columns",
            "multischolar_invalid_prot_dia_da_table"
        )
    }
    keys <- table[c(protein_id_column, "comparison")]
    if (any(!stats::complete.cases(keys)) || anyDuplicated(keys) > 0L) {
        protDiaDaArtifactAbort(
            "DIA-NN DA long table has missing or duplicate result identities",
            "multischolar_invalid_prot_dia_da_identity"
        )
    }
    observed <- unique(as.character(table$comparison))
    expected <- c(metadata$friendly_name, metadata$full_format)
    if (length(observed) > 1L || !all(observed %in% expected)) {
        protDiaDaArtifactAbort(
            "DIA-NN DA long table comparison differs from its contrast",
            "multischolar_invalid_prot_dia_da_identity"
        )
    }
    invisible(table)
}

protDiaDaQuerySpecification <- function(ref, table, protein_id_column) {
    projections <- names(table)
    filters <- list()
    if (is.character(table[[protein_id_column]])) {
        filters$protein_search <- list(
            column = protein_id_column,
            type = "character",
            operators = "contains"
        )
    }
    filters$q_value <- list(
        column = "fdr_qvalue",
        type = "double",
        operators = c("lt", "lte", "gt", "gte", "between", "is_missing")
    )
    filters$log2fc_min <- list(
        column = "log2FC",
        type = "double",
        operators = c("gt", "gte")
    )
    filters$log2fc_max <- list(
        column = "log2FC",
        type = "double",
        operators = c("lt", "lte")
    )
    sorts <- list(
        row_order = list(
            column = "comparison",
            transform = "identity",
            directions = "asc"
        ),
        protein = list(
            column = protein_id_column,
            transform = "identity",
            directions = c("asc", "desc")
        ),
        q_value = list(
            column = "fdr_qvalue",
            transform = "identity",
            directions = c("asc", "desc")
        ),
        effect = list(
            column = "log2FC",
            transform = "identity",
            directions = c("asc", "desc")
        ),
        absolute_effect = list(
            column = "log2FC",
            transform = "absolute",
            directions = c("asc", "desc")
        )
    )
    newArtifactQueryPageSpecification(
        query_id = paste0(
            "proteomics.",
            ref$logical_key$workflow_slug,
            ".da.",
            ref$artifact_id,
            ".v1"
        ),
        state_role = ref$logical_key$state_role,
        projections = projections,
        filters = filters,
        sorts = sorts,
        default_sort = list(sort_id = "row_order", direction = "asc"),
        max_rows = .PROT_DIA_DA_MAX_QUERY_ROWS,
        max_bytes = .PROT_DIA_DA_MAX_QUERY_BYTES
    )
}

protDiaDaTableRole <- function(contrast_id, table_name) {
    safe_name <- gsub("[^A-Za-z0-9_]", "_", table_name)
    paste0("da_", sub("^contrast_", "", contrast_id), "_", safe_name)
}

protDiaDaEncodeTable <- function(
    table,
    table_name,
    contrast_id,
    protein_id_column
) {
    stable_key <- if (identical(table_name, "da_proteins_long")) {
        c(protein_id_column, "comparison")
    } else {
        NULL
    }
    encodeArtifactTable(
        table,
        stable_key = stable_key,
        owner = paste0("DIA-NN DA ", contrast_id, " ", table_name)
    )
}

protDiaDaWriteTable <- function(
    store,
    run_id,
    contrast_id,
    table_name,
    table,
    protein_id_column
) {
    encoded <- protDiaDaEncodeTable(
        table,
        table_name,
        contrast_id,
        protein_id_column
    )
    identity <- store$labels
    logical_key <- list(
        project_id = store$project_id,
        omic_type = identity$omic_type,
        workflow_slug = identity$workflow_slug,
        stage_id = "differential_abundance",
        state_role = protDiaDaTableRole(contrast_id, table_name),
        generation_id = run_id
    )
    ref <- artifactStoreWriteParquet(
        store,
        encoded,
        logical_key,
        provenance_ids = run_id
    )
    list(
        table_name = table_name,
        ref = unclass(ref),
        semantic_digest = encoded$metadata$semantic_digest,
        rows = as.integer(nrow(table)),
        columns = as.integer(ncol(table))
    )
}

protDiaDaModelDigest <- function(model) {
    digest::digest(
        model,
        algo = .artifactHashAlgorithm,
        serialize = TRUE,
        ascii = FALSE,
        serializeVersion = 3L
    )
}

protDiaDaWriteModel <- function(
    model,
    project_root,
    relative_path,
    save_rds_fn = saveRDS
) {
    path <- artifactResolveContainedPath(project_root, relative_path)
    if (file.exists(path) || dir.exists(path)) {
        protDiaDaArtifactAbort(
            "immutable DIA-NN DA model already exists",
            "multischolar_prot_dia_da_model_exists"
        )
    }
    parent <- dirname(path)
    if (!dir.exists(parent) &&
        !dir.create(parent, recursive = TRUE, showWarnings = FALSE)) {
        protDiaDaArtifactAbort(
            "DIA-NN DA model directory could not be created",
            "multischolar_prot_dia_da_write_failed"
        )
    }
    temporary <- file.path(
        parent,
        paste0(".", basename(path), ".", artifactOpaqueId("tmp"), ".tmp")
    )
    on.exit(if (file.exists(temporary)) unlink(temporary, force = FALSE), add = TRUE)
    save_rds_fn(model, temporary, version = 3)
    restored <- readRDS(temporary)
    semantic_digest <- protDiaDaModelDigest(model)
    if (!identical(protDiaDaModelDigest(restored), semantic_digest)) {
        protDiaDaArtifactAbort(
            "DIA-NN DA model changed during RDS serialization",
            "multischolar_inexact_prot_dia_da_model"
        )
    }
    if (!isTRUE(file.rename(temporary, path))) {
        protDiaDaArtifactAbort(
            "DIA-NN DA model could not be atomically published",
            "multischolar_prot_dia_da_publish_failed"
        )
    }
    list(
        schema = .PROT_DIA_DA_MODEL_SCHEMA,
        schema_version = .PROT_DIA_DA_MODEL_VERSION,
        relative_path = relative_path,
        byte_digest = artifactByteDigest(path),
        semantic_digest = semantic_digest,
        class = unname(class(model))
    )
}

protDiaDaValidateModel <- function(model_ref, project_root) {
    required <- c(
        "schema", "schema_version", "relative_path", "byte_digest",
        "semantic_digest", "class"
    )
    valid <- is.list(model_ref) && identical(names(model_ref), required) &&
        identical(model_ref$schema, .PROT_DIA_DA_MODEL_SCHEMA) &&
        identical(
            workflowStateVersionValue(model_ref$schema_version),
            .PROT_DIA_DA_MODEL_VERSION
        ) && (is.character(model_ref$class) || is.list(model_ref$class)) &&
        length(model_ref$class) > 0L
    if (!isTRUE(valid)) {
        protDiaDaArtifactAbort(
            "DIA-NN DA model reference is malformed",
            "multischolar_invalid_prot_dia_da_model"
        )
    }
    model_ref$class <- unlist(model_ref$class, use.names = FALSE)
    artifactRefValidateDigest(model_ref$byte_digest, "model byte digest")
    artifactRefValidateDigest(
        model_ref$semantic_digest, "model semantic digest"
    )
    path <- artifactResolveContainedPath(
        project_root, model_ref$relative_path, must_exist = TRUE
    )
    if (!identical(artifactByteDigest(path), model_ref$byte_digest) ||
        !identical(protDiaDaModelDigest(readRDS(path)),
            model_ref$semantic_digest)) {
        protDiaDaArtifactAbort(
            "DIA-NN DA model fingerprint validation failed",
            "multischolar_prot_dia_da_model_digest_mismatch"
        )
    }
    model_ref$schema_version <- .PROT_DIA_DA_MODEL_VERSION
    model_ref
}

protDiaDaSummary <- function(table, parameters) {
    q_value <- table$fdr_qvalue
    effect <- table$log2FC
    significant <- q_value < parameters$da_q_val_thresh
    list(
        rows = as.integer(nrow(table)),
        significant = as.integer(sum(significant, na.rm = TRUE)),
        up = as.integer(sum(
            significant & effect > parameters$treat_lfc_cutoff,
            na.rm = TRUE
        )),
        down = as.integer(sum(
            significant & effect < -parameters$treat_lfc_cutoff,
            na.rm = TRUE
        )),
        q_value_threshold = parameters$da_q_val_thresh,
        fold_change_cutoff = parameters$treat_lfc_cutoff
    )
}

protDiaDaValidateTableEntry <- function(entry, store) {
    required <- c(
        "table_name", "ref", "semantic_digest", "rows", "columns"
    )
    if (!is.list(entry) || !identical(names(entry), required) ||
        !workflowCapabilityScalarString(entry$table_name)) {
        protDiaDaArtifactAbort(
            "DIA-NN DA table entry is malformed",
            "multischolar_invalid_prot_dia_da_manifest"
        )
    }
    ref <- artifactStoreNormalizeRef(entry$ref)
    sidecar <- artifactStoreReadSidecar(
        store,
        artifactStoreManagedPaths(
            store, ref$logical_key, ref$artifact_id
        )$sidecar,
        validate_payload = TRUE
    )
    if (!identical(sidecar$artifact_ref, ref) ||
        !identical(ref$hash_policy$semantic$digest, entry$semantic_digest) ||
        !identical(as.integer(entry$rows), ref$shape$rows) ||
        !identical(as.integer(entry$columns), ref$shape$columns)) {
        protDiaDaArtifactAbort(
            "DIA-NN DA table entry differs from its immutable artifact",
            "multischolar_prot_dia_da_table_mismatch"
        )
    }
    entry$ref <- unclass(ref)
    entry$rows <- as.integer(entry$rows)
    entry$columns <- as.integer(entry$columns)
    entry
}

protDiaDaValidateContrastManifest <- function(manifest, store = NULL) {
    required <- c(
        "schema", "schema_version", "project_id", "workflow_id", "run_id",
        "source_generation_id", "contrast", "tables", "model", "summary",
        "query_specification", "manifest_digest"
    )
    valid <- is.list(manifest) && identical(names(manifest), required) &&
        identical(manifest$schema, .PROT_DIA_DA_CONTRAST_SCHEMA) &&
        identical(
            workflowStateVersionValue(manifest$schema_version),
            .PROT_DIA_DA_CONTRAST_VERSION
        ) && all(vapply(
            manifest[c(
                "project_id", "workflow_id", "run_id", "source_generation_id"
            )],
            workflowCapabilityScalarString,
            logical(1)
        )) && is.list(manifest$contrast) &&
        identical(
            names(manifest$contrast),
            c("contrast_id", "contrast_name", "full_format", "friendly_name")
        ) && all(vapply(
            manifest$contrast,
            workflowCapabilityScalarString,
            logical(1)
        )) && is.list(manifest$tables) &&
        "da_proteins_long" %in% names(manifest$tables)
    if (!isTRUE(valid)) {
        protDiaDaArtifactAbort(
            "DIA-NN DA contrast manifest is malformed or unsupported",
            "multischolar_invalid_prot_dia_da_manifest"
        )
    }
    artifactWorkflowStateAssertSafeMetadata(
        manifest, "DIA-NN DA contrast manifest"
    )
    expected_digest <- protDiaDaJsonDigest(manifest)
    if (!identical(manifest$manifest_digest, expected_digest)) {
        protDiaDaArtifactAbort(
            "DIA-NN DA contrast manifest digest differs",
            "multischolar_prot_dia_da_manifest_digest_mismatch"
        )
    }
    manifest$query_specification <- validateArtifactQueryPageSpecification(
        manifest$query_specification
    )
    manifest$summary <- protDiaDaValidateSummary(manifest$summary)
    if (!is.null(store)) {
        manifest$tables <- lapply(
            manifest$tables, protDiaDaValidateTableEntry, store = store
        )
        table_names <- vapply(
            manifest$tables, `[[`, character(1), "table_name"
        )
        if (!identical(unname(table_names), names(manifest$tables))) {
            protDiaDaArtifactAbort(
                "DIA-NN DA table names disagree with their manifest keys",
                "multischolar_prot_dia_da_table_mismatch"
            )
        }
        manifest$model <- protDiaDaValidateModel(
            manifest$model, store$project_root
        )
        long_ref <- manifest$tables$da_proteins_long$ref
        if (!identical(
            manifest$query_specification$state_role,
            long_ref$logical_key$state_role
        )) {
            protDiaDaArtifactAbort(
                "DIA-NN DA query specification targets the wrong table",
                "multischolar_prot_dia_da_query_mismatch"
            )
        }
    }
    manifest$schema_version <- .PROT_DIA_DA_CONTRAST_VERSION
    manifest
}

protDiaDaPersistContrast <- function(
    workflow_data,
    run_id,
    source,
    metadata,
    result,
    parameters,
    failure_injector = NULL
) {
    context <- workflow_data$workflow_context
    identity <- context$getIdentity()
    store <- protDiaDaStore(context)
    contrast_id <- protDiaDaContrastId(metadata)
    paths <- protDiaDaPaths(context, run_id, contrast_id)
    manifest_path <- artifactResolveContainedPath(
        context$getProjectRoot(), paths$contrast_manifest
    )
    if (file.exists(manifest_path)) {
        manifest <- protDiaDaReadJson(
            manifest_path,
            function(value) protDiaDaValidateContrastManifest(value, store)
        )
        result_tables <- protDiaDaResultTables(result)
        protein_id_column <- result$theObject@protein_id_column
        expected <- vapply(names(result_tables), function(table_name) {
            protDiaDaEncodeTable(
                result_tables[[table_name]],
                table_name,
                contrast_id,
                protein_id_column
            )$metadata$semantic_digest
        }, character(1))
        observed <- vapply(
            manifest$tables,
            `[[`,
            character(1),
            "semantic_digest"
        )
        if (!identical(names(observed), names(expected)) ||
            !identical(unname(observed), unname(expected)) ||
            !identical(
                manifest$model$semantic_digest,
                protDiaDaModelDigest(result$contrasts_results)
            )) {
            protDiaDaArtifactAbort(
                "resumed DIA-NN DA contrast differs from committed content",
                "multischolar_prot_dia_da_resume_mismatch"
            )
        }
        return(manifest)
    }
    tables <- protDiaDaResultTables(result)
    if (!"da_proteins_long" %in% names(tables)) {
        protDiaDaArtifactAbort(
            "DIA-NN DA contrast has no long result table",
            "multischolar_invalid_prot_dia_da_table"
        )
    }
    protein_id_column <- result$theObject@protein_id_column
    protDiaDaValidateLongTable(
        tables$da_proteins_long, protein_id_column, metadata
    )
    table_entries <- lapply(names(tables), function(table_name) {
        protDiaDaWriteTable(
            store,
            run_id,
            contrast_id,
            table_name,
            tables[[table_name]],
            protein_id_column
        )
    })
    names(table_entries) <- names(tables)
    long_ref <- artifactStoreNormalizeRef(
        table_entries$da_proteins_long$ref
    )
    query_specification <- protDiaDaQuerySpecification(
        long_ref, tables$da_proteins_long, protein_id_column
    )
    model <- protDiaDaWriteModel(
        result$contrasts_results,
        context$getProjectRoot(),
        paths$model
    )
    contrast <- c(list(contrast_id = contrast_id), metadata)
    manifest <- list(
        schema = .PROT_DIA_DA_CONTRAST_SCHEMA,
        schema_version = .PROT_DIA_DA_CONTRAST_VERSION,
        project_id = identity$project_id,
        workflow_id = identity$workflow_id,
        run_id = run_id,
        source_generation_id = source$generation_id,
        contrast = contrast,
        tables = table_entries,
        model = model,
        summary = protDiaDaSummary(tables$da_proteins_long, parameters),
        query_specification = query_specification,
        manifest_digest = NULL
    )
    manifest$manifest_digest <- protDiaDaJsonDigest(manifest)
    protDiaDaWriteJson(
        manifest,
        manifest_path,
        function(value) protDiaDaValidateContrastManifest(value, store)
    )
    artifactStoreInvokeFailure(
        failure_injector,
        "after_contrast_commit",
        list(run_id = run_id, contrast_id = contrast_id)
    )
    protDiaDaValidateContrastManifest(manifest, store)
}

prepareProtDiaDaArtifactRun <- function(
    workflow_data,
    contrasts,
    results,
    parameters,
    failure_injector = NULL,
    now = Sys.time(),
    eligibility_fn = protDiaDaArtifactEligible
) {
    if (!eligibility_fn(workflow_data)) {
        return(list(enabled = FALSE, index = NULL, manifest = NULL))
    }
    contrasts <- normaliseProtDaContrastsTable(contrasts)
    if (!identical(names(results), contrasts$contrasts)) {
        protDiaDaArtifactAbort(
            "DIA-NN DA result order differs from requested contrasts",
            "multischolar_prot_dia_da_contrast_mismatch"
        )
    }
    source <- protDiaDaSourceBinding(workflow_data)
    software <- protDiaDaSoftware()
    run_id <- protDiaDaRunId(
        source, protDiaDaRunContrasts(contrasts), parameters, software
    )
    context <- workflow_data$workflow_context
    paths <- protDiaDaPaths(context, run_id)
    run_path <- artifactResolveContainedPath(
        context$getProjectRoot(), paths$run_manifest
    )
    if (file.exists(run_path)) {
        manifest <- protDiaDaReadJson(
            run_path,
            function(value) protDiaDaValidateRunManifest(value, context)
        )
        return(list(
            enabled = TRUE,
            resumed = TRUE,
            index = protDiaDaRunIndex(
                manifest, paths$run_manifest, context
            ),
            manifest = manifest,
            manifest_path = run_path
        ))
    }
    contrast_manifests <- lapply(seq_len(nrow(contrasts)), function(index) {
        metadata <- protDiaDaContrastMetadata(contrasts, index)
        protDiaDaPersistContrast(
            workflow_data,
            run_id,
            source,
            metadata,
            results[[index]],
            parameters,
            failure_injector
        )
    })
    entries <- lapply(contrast_manifests, function(contrast) {
        contrast_paths <- protDiaDaPaths(
            context, run_id, contrast$contrast$contrast_id
        )
        list(
            contrast_id = contrast$contrast$contrast_id,
            contrast_name = contrast$contrast$contrast_name,
            full_format = contrast$contrast$full_format,
            friendly_name = contrast$contrast$friendly_name,
            manifest_relative_path = contrast_paths$contrast_manifest,
            manifest_digest = contrast$manifest_digest
        )
    })
    identity <- context$getIdentity()
    descriptor <- findArtifactWorkflowDescriptor(
        identity, artifactWorkflowDescriptorCatalogue()
    )
    manifest <- list(
        schema = .PROT_DIA_DA_RUN_SCHEMA,
        schema_version = .PROT_DIA_DA_RUN_VERSION,
        project_id = identity$project_id,
        workflow_id = identity$workflow_id,
        descriptor_contract = list(
            descriptor_id = descriptor$descriptor_id,
            descriptor_version = descriptor$descriptor_version,
            descriptor_digest = descriptor$descriptor_digest
        ),
        run_id = run_id,
        source = source,
        parameters = parameters,
        parameters_digest = artifactSemanticDigest(parameters),
        software = software,
        contrasts = entries,
        created_at = format(now, "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC"),
        manifest_digest = NULL
    )
    manifest$manifest_digest <- protDiaDaJsonDigest(manifest)
    protDiaDaWriteJson(
        manifest,
        run_path,
        function(value) protDiaDaValidateRunManifest(value, context)
    )
    manifest <- protDiaDaValidateRunManifest(manifest, context)
    list(
        enabled = TRUE,
        resumed = FALSE,
        index = protDiaDaRunIndex(manifest, paths$run_manifest, context),
        manifest = manifest,
        manifest_path = run_path
    )
}

publishProtDiaDaArtifactRun <- function(
    workflow_data,
    prepared,
    failure_injector = NULL
) {
    if (!isTRUE(prepared$enabled)) return(prepared)
    context <- workflow_data$workflow_context
    current_generation <- workflow_data$state_manager$getCurrentGenerationId()
    if (!identical(
        current_generation,
        prepared$manifest$source$generation_id
    )) {
        protDiaDaArtifactAbort(
            "DIA-NN DA source generation changed before publication",
            "multischolar_prot_dia_da_source_mismatch"
        )
    }
    current_path <- artifactResolveContainedPath(
        context$getProjectRoot(), protDiaDaPaths(context)$current
    )
    protDiaDaWriteJson(
        prepared$manifest,
        current_path,
        function(value) protDiaDaValidateRunManifest(value, context),
        replace = TRUE,
        failure_injector = failure_injector,
        failure_stage = "before_da_current_publication"
    )
    prepared$current_path <- current_path
    prepared$published <- TRUE
    prepared
}
