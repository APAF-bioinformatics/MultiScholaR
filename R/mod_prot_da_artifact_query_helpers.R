# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

protDiaDaResolveContrasts <- function(workflow_data, da_data = NULL) {
    if (protDiaDaArtifactWorkflow(workflow_data)) {
        contrasts <- workflow_data$contrasts_tbl
        if (is.null(contrasts) || !is.data.frame(contrasts) ||
            nrow(contrasts) == 0L) {
            protDiaDaArtifactAbort(
                "artifact-mode DIA-NN DA has no explicit workflow contrasts",
                "multischolar_missing_prot_dia_da_contrasts"
            )
        }
        return(normaliseProtDaContrastsTable(contrasts))
    }
    if (exists("contrasts_tbl", envir = .GlobalEnv, inherits = FALSE)) {
        return(normaliseProtDaContrastsTable(get(
            "contrasts_tbl", envir = .GlobalEnv, inherits = FALSE
        )))
    }
    available <- if (is.null(da_data)) NULL else da_data$contrasts_available
    normaliseProtDaContrastsTable(data.frame(
        contrasts = available,
        stringsAsFactors = FALSE
    ))
}

protDiaDaResolveAnnotations <- function(workflow_data) {
    if (protDiaDaArtifactWorkflow(workflow_data)) {
        return(workflow_data$uniprot_dat_cln)
    }
    if (exists("uniprot_dat_cln", envir = .GlobalEnv, inherits = FALSE)) {
        return(get("uniprot_dat_cln", envir = .GlobalEnv, inherits = FALSE))
    }
    NULL
}

protDiaDaContrastsVector <- function(workflow_data, da_data = NULL) {
    contrasts <- protDiaDaResolveContrasts(workflow_data, da_data)
    if ("comparison" %in% names(contrasts)) return(contrasts$comparison)
    if ("contrasts" %in% names(contrasts)) return(contrasts$contrasts)
    contrasts[[1L]]
}


protDiaDaValidateRunSource <- function(source, store) {
    required <- c(
        "generation_id", "state_name", "manifest_relative_path",
        "manifest_digest", "state_semantic_digest", "stage_refs"
    )
    if (!is.list(source) || !identical(names(source), required)) {
        protDiaDaArtifactAbort(
            "DIA-NN DA source binding is malformed",
            "multischolar_prot_dia_da_source_mismatch"
        )
    }
    manifest <- artifactWorkflowStateReadManifest(
        store, source$manifest_relative_path
    )
    valid <- identical(manifest$generation_id, source$generation_id) &&
        identical(manifest$logical_name, source$state_name) &&
        identical(manifest$manifest_digest, source$manifest_digest) &&
        identical(
            manifest$data$semantic_digest,
            source$state_semantic_digest
        )
    if (!isTRUE(valid)) {
        protDiaDaArtifactAbort(
            "DIA-NN DA quantitative source fingerprint differs",
            "multischolar_prot_dia_da_source_mismatch"
        )
    }
    protDiaSessionValidateStage(
        store,
        source$stage_refs$import,
        .PROT_DIA_SESSION_IMPORT_ROLES,
        "import"
    )
    protDiaSessionValidateStage(
        store,
        source$stage_refs$design,
        .PROT_DIA_SESSION_DESIGN_ROLES,
        "design"
    )
    invisible(TRUE)
}

protDiaDaRunIndex <- function(
    manifest,
    manifest_relative_path,
    context = NULL
) {
    entries <- lapply(manifest$contrasts, function(entry) {
        contrast <- entry$manifest
        if (is.null(contrast)) {
            if (is.null(context)) {
                protDiaDaArtifactAbort(
                    "DIA-NN DA index requires a verified project context",
                    "multischolar_invalid_prot_dia_da_index"
                )
            }
            store <- protDiaDaStore(context)
            path <- artifactResolveContainedPath(
                context$getProjectRoot(),
                entry$manifest_relative_path,
                must_exist = TRUE
            )
            contrast <- protDiaDaReadJson(
                path,
                function(value) {
                    protDiaDaValidateContrastManifest(value, store)
                }
            )
        }
        list(
            contrast_id = entry$contrast_id,
            contrast_name = entry$contrast_name,
            full_format = entry$full_format,
            friendly_name = entry$friendly_name,
            manifest_relative_path = entry$manifest_relative_path,
            manifest_digest = entry$manifest_digest,
            long_table = contrast$tables$da_proteins_long,
            query_specification = contrast$query_specification,
            summary = contrast$summary
        )
    })
    structure(
        list(
            schema = .PROT_DIA_DA_INDEX_SCHEMA,
            schema_version = .PROT_DIA_DA_INDEX_VERSION,
            backend = "artifact",
            run_id = manifest$run_id,
            source_generation_id = manifest$source$generation_id,
            manifest_relative_path = manifest_relative_path,
            manifest_digest = manifest$manifest_digest,
            contrasts = entries
        ),
        class = c("MultiScholaRProtDiaDaIndex", "list")
    )
}

isProtDiaDaArtifactIndex <- function(value) {
    inherits(value, "MultiScholaRProtDiaDaIndex") &&
        identical(value$schema, .PROT_DIA_DA_INDEX_SCHEMA) &&
        identical(value$schema_version, .PROT_DIA_DA_INDEX_VERSION) &&
        identical(value$backend, "artifact")
}

restoreProtDiaDaArtifactIndex <- function(workflow_data) {
    if (!protDiaDaArtifactEligible(workflow_data, "queries")) return(NULL)
    context <- workflow_data$workflow_context
    current_path <- artifactResolveContainedPath(
        context$getProjectRoot(), protDiaDaPaths(context)$current
    )
    if (!file.exists(current_path)) return(NULL)
    manifest <- protDiaDaReadJson(
        current_path,
        function(value) protDiaDaValidateRunManifest(value, context)
    )
    if (!identical(
        manifest$source$generation_id,
        workflow_data$state_manager$getCurrentGenerationId()
    )) {
        return(NULL)
    }
    paths <- protDiaDaPaths(context, manifest$run_id)
    protDiaDaRunIndex(manifest, paths$run_manifest, context)
}

applyProtDiaDaArtifactIndex <- function(workflow_data, da_data) {
    index <- restoreProtDiaDaArtifactIndex(workflow_data)
    if (is.null(index)) return(FALSE)
    da_data$da_results_list <- index
    da_data$analysis_complete <- TRUE
    workflow_data$da_analysis_results_list <- index
    status <- workflow_data$tab_status
    status$differential_expression <- "complete"
    status$differential_abundance <- "complete"
    status$enrichment_analysis <- "pending"
    workflow_data$tab_status <- status
    TRUE
}

protDiaDaIndexEntry <- function(index, contrast) {
    if (!isProtDiaDaArtifactIndex(index) ||
        !workflowCapabilityScalarString(contrast)) {
        protDiaDaArtifactAbort(
            "DIA-NN DA artifact index request is invalid",
            "multischolar_invalid_prot_dia_da_index"
        )
    }
    matches <- vapply(index$contrasts, function(entry) {
        contrast %in% c(
            entry$contrast_name, entry$full_format, entry$friendly_name
        )
    }, logical(1))
    if (sum(matches) != 1L) {
        protDiaDaArtifactAbort(
            "DIA-NN DA contrast selection is absent or ambiguous",
            "multischolar_unknown_prot_dia_da_contrast"
        )
    }
    index$contrasts[[which(matches)]]
}

queryProtDiaDaPage <- function(
    workflow_data,
    index,
    contrast,
    projections = NULL,
    filters = list(),
    sort_id = NULL,
    direction = NULL,
    cursor = NULL,
    limit = 100L,
    resource_policy = NULL
) {
    if (!protDiaDaArtifactEligible(workflow_data, "queries")) {
        protDiaDaArtifactAbort(
            "DIA-NN DA server-side queries are disabled or ineligible",
            "multischolar_prot_dia_da_queries_disabled"
        )
    }
    entry <- protDiaDaIndexEntry(index, contrast)
    store <- protDiaDaStore(workflow_data$workflow_context)
    queryArtifactPage(
        store = store,
        ref = entry$long_table$ref,
        specification = entry$query_specification,
        projections = projections,
        filters = filters,
        sort_id = sort_id,
        direction = direction,
        cursor = cursor,
        limit = limit,
        resource_policy = resource_policy
    )
}

protDiaDaCompleteSelectedTable <- function(
    workflow_data,
    index,
    contrast,
    projections = NULL,
    filters = list(),
    sort_id = NULL,
    direction = NULL,
    limit = NULL,
    resource_policy = NULL
) {
    entry <- protDiaDaIndexEntry(index, contrast)
    rows <- as.integer(entry$long_table$rows)
    if (is.null(limit)) limit <- rows
    if (limit < rows && length(filters) == 0L) {
        protDiaDaArtifactAbort(
            "selected DIA-NN DA contrast exceeds the complete-view row bound",
            "multischolar_prot_dia_da_complete_view_limit"
        )
    }
    page <- queryProtDiaDaPage(
        workflow_data,
        index,
        contrast,
        projections,
        filters,
        sort_id,
        direction,
        cursor = NULL,
        limit = limit,
        resource_policy = resource_policy
    )
    if (isTRUE(page$has_more)) {
        protDiaDaArtifactAbort(
            "selected DIA-NN DA view would truncate scientific results",
            "multischolar_prot_dia_da_complete_view_limit"
        )
    }
    page$data
}

protDiaDaSelectedResults <- function(
    workflow_data,
    da_data,
    contrast,
    view = c("volcano", "heatmap"),
    q_value_threshold = NULL,
    top_n = NULL
) {
    view <- match.arg(view)
    index <- da_data$da_results_list
    entry <- protDiaDaIndexEntry(index, contrast)
    protein_id <- da_data$current_s4_object@protein_id_column
    if (identical(view, "heatmap")) {
        projections <- intersect(
            c(protein_id, "comparison", "log2FC", "fdr_qvalue", "gene_name"),
            entry$query_specification$projections
        )
        filters <- list(q_value = list(
            operator = "lt", value = as.double(q_value_threshold)
        ))
        page <- queryProtDiaDaPage(
            workflow_data,
            index,
            contrast,
            projections = projections,
            filters = filters,
            sort_id = "absolute_effect",
            direction = "desc",
            limit = as.integer(top_n)
        )
        table <- page$data
    } else {
        table <- protDiaDaCompleteSelectedTable(
            workflow_data,
            index,
            contrast
        )
    }
    list(
        da_proteins_long = table,
        theObject = da_data$current_s4_object
    )
}

protDiaDaTableFilters <- function(significance, q_value, lfc) {
    filters <- list()
    if (significance %in% c("significant", "up", "down")) {
        filters$q_value <- list(operator = "lt", value = as.double(q_value))
    }
    if (identical(significance, "up")) {
        filters$log2fc_min <- list(operator = "gt", value = as.double(lfc))
    }
    if (identical(significance, "down")) {
        filters$log2fc_max <- list(operator = "lt", value = -as.double(lfc))
    }
    filters
}

protDiaDaArtifactSummary <- function(
    workflow_data,
    index,
    contrast,
    q_value,
    lfc
) {
    entry <- protDiaDaIndexEntry(index, contrast)
    summary <- entry$summary
    if (identical(as.double(summary$q_value_threshold), as.double(q_value)) &&
        identical(as.double(summary$fold_change_cutoff), as.double(lfc))) {
        return(summary)
    }
    table <- protDiaDaCompleteSelectedTable(
        workflow_data,
        index,
        contrast,
        projections = c("fdr_qvalue", "log2FC")
    )
    protDiaDaSummary(
        table,
        list(
            da_q_val_thresh = as.double(q_value),
            treat_lfc_cutoff = as.double(lfc)
        )
    )
}
