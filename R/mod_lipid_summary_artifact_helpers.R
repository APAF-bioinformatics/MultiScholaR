# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

lipidSummaryArtifactMode <- function() {
    match.arg(
        getOption("multischolar.lipidomics.summary_artifacts", "enabled"),
        c("enabled", "disabled")
    )
}

lipidSummaryArtifactEligible <- function(workflow_data) {
    !identical(lipidSummaryArtifactMode(), "disabled") &&
        lipidQcWorkflowData(workflow_data) &&
        inherits(workflow_data$state_manager, "ArtifactWorkflowState") &&
        inherits(workflow_data$workflow_context, "WorkflowContext") &&
        workflow_data$workflow_context$isBound() &&
        identical(
            workflow_data$workflow_context$getStorageDecision()$effective_backend,
            "artifact"
        )
}

prepareLipidSummaryDependencies <- function(
    workflow_data,
    project_dirs,
    omic_type
) {
    if (!lipidSummaryArtifactEligible(workflow_data)) {
        return(list(enabled = FALSE, reason = "artifact_summary_disabled"))
    }
    object <- workflow_data$state_manager$getState()
    if (!methods::is(object, "LipidomicsAssayData") ||
        !identical(methods::validObject(object, test = TRUE), TRUE)) {
        workflowSessionAbort(
            "lipidomics summary requires a valid final assay object",
            "multischolar_missing_lipid_summary_dependency"
        )
    }
    da <- workflow_data$artifact_stage_results$differential_abundance
    if (!is.list(da) || !isTRUE(da$committed)) {
        workflowSessionAbort(
            "lipidomics summary requires current DA artifacts",
            "multischolar_missing_lipid_summary_dependency"
        )
    }
    if (!is.list(project_dirs) || !omic_type %in% names(project_dirs)) {
        workflowSessionAbort(
            "lipidomics summary requires explicit project directories",
            "multischolar_missing_lipid_summary_dependency"
        )
    }
    paths <- project_dirs[[omic_type]]
    context <- workflow_data$workflow_context
    required_paths <- c("base_dir", "source_dir", "results_summary_dir")
    if (!all(vapply(required_paths, function(name) {
        workflowCapabilityScalarString(paths[[name]]) &&
            artifactPathIsContained(
                normalizePath(paths[[name]], winslash = "/", mustWork = FALSE),
                context$getProjectRoot()
            )
    }, logical(1)))) {
        workflowSessionAbort(
            "lipidomics summary paths are incomplete or outside the project",
            "multischolar_lipid_summary_path_mismatch"
        )
    }
    contrasts <- workflow_data$contrasts_tbl
    if (!is.data.frame(contrasts) || nrow(contrasts) == 0L) {
        workflowSessionAbort(
            "lipidomics summary requires contrasts",
            "multischolar_missing_lipid_summary_dependency"
        )
    }
    dependencies <- new.env(parent = emptyenv())
    dependencies$enabled <- TRUE
    dependencies$object <- object
    dependencies$da <- da
    dependencies$contrasts <- contrasts
    dependencies$config <- workflow_data$config_list
    dependencies$paths <- paths
    dependencies$state_generation_id <-
        workflow_data$state_manager$getCurrentGenerationId()
    dependencies$assay_order <- names(object@lipid_data)
    dependencies$evidence <- list(
        state_generation_id = dependencies$state_generation_id,
        state_digest = lipidNormSafeDigest(object),
        da_run_id = da$run_id,
        da_manifest_digest = da$manifest$manifest_digest,
        assay_order = dependencies$assay_order,
        design_digest = artifactSemanticDigest(object@design_matrix),
        contrasts_digest = artifactSemanticDigest(contrasts),
        args_digest = lipidNormSafeDigest(object@args)
    )
    dependencies
}

releaseLipidSummaryDependencies <- function(dependencies) {
    if (is.environment(dependencies)) {
        rm(list = ls(dependencies, all.names = TRUE), envir = dependencies)
    }
    invisible(TRUE)
}

recordLipidSummaryProduct <- function(dependencies, path, product_type) {
    if (!is.environment(dependencies) || !file.exists(path)) {
        return(invisible(NULL))
    }
    evidence <- dependencies$evidence
    record <- list(
        schema = "multischolar.lipidomics_summary_product",
        schema_version = 1L,
        product_type = product_type,
        relative_path = workflowSessionProjectRelativePath(
            path,
            dependencies$paths$base_dir
        ),
        byte_digest = artifactByteDigest(path),
        bytes = unname(as.numeric(file.info(path)$size)),
        dependencies = evidence,
        created_at = artifactRefUtcNow(),
        manifest_digest = NULL
    )
    record$manifest_digest <- workflowSessionContentDigest(record)
    metadata_dir <- file.path(
        dependencies$paths$results_summary_dir,
        "artifact_metadata"
    )
    dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)
    metadata_path <- file.path(
        metadata_dir,
        paste0(product_type, ".json")
    )
    temporary <- paste0(metadata_path, ".", artifactOpaqueId("tmp"), ".tmp")
    on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
    writeLines(workflowSessionJson(record), temporary)
    if (!isTRUE(file.rename(temporary, metadata_path))) {
        workflowSessionAbort(
            "lipidomics summary product metadata publication failed",
            "multischolar_lipid_summary_publish_failed"
        )
    }
    invisible(metadata_path)
}

lipidSummaryGlobalOwnershipAllowed <- function(workflow_data) {
    !lipidSummaryArtifactEligible(workflow_data)
}
