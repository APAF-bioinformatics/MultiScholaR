# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

#' Resolve the independent non-DIA DA artifact switch
#' @param kind Persistence or server-side query direction.
#' @return One of `"enabled"` or `"disabled"`.
#' @noRd
protNonDiaDaArtifactMode <- function(kind = c("persistence", "queries")) {
    kind <- match.arg(kind)
    option <- paste0("multischolar.prot_nondia.da_", kind)
    match.arg(getOption(option, "enabled"), c("enabled", "disabled"))
}

#' Test whether workflow state is an exact non-DIA artifact DA workflow
#' @param workflow_data Mutable proteomics workflow state.
#' @return A scalar logical.
#' @noRd
protNonDiaDaArtifactWorkflow <- function(workflow_data) {
    descriptor <- protNonDiaNormDescriptor(workflow_data)
    manager <- tryCatch(workflow_data$state_manager, error = \(...) NULL)
    !is.null(descriptor) && inherits(manager, "ArtifactWorkflowState") &&
        identical(
            workflowStateType(manager),
            descriptor$identity$workflow_type
        ) && artifactStageCoordinatorOwned(workflow_data, descriptor)
}

#' Test exact non-DIA DA artifact eligibility
#' @param workflow_data Mutable proteomics workflow state.
#' @param kind Persistence or server-side query direction.
#' @return A scalar logical.
#' @noRd
protNonDiaDaArtifactEligible <- function(
    workflow_data,
    kind = "persistence"
) {
    identical(protNonDiaDaArtifactMode(kind), "enabled") &&
        protNonDiaDaArtifactWorkflow(workflow_data)
}

#' Test any exact supported proteomics artifact DA workflow
#' @param workflow_data Mutable proteomics workflow state.
#' @return A scalar logical.
#' @noRd
protDaArtifactWorkflow <- function(workflow_data) {
    protDiaDaArtifactWorkflow(workflow_data) ||
        protNonDiaDaArtifactWorkflow(workflow_data)
}

#' Test descriptor-dispatched DA artifact eligibility
#' @param workflow_data Mutable proteomics workflow state.
#' @param kind Persistence or server-side query direction.
#' @return A scalar logical.
#' @noRd
protDaArtifactEligible <- function(workflow_data, kind = "persistence") {
    if (protDiaDaArtifactWorkflow(workflow_data)) {
        return(protDiaDaArtifactEligible(workflow_data, kind))
    }
    protNonDiaDaArtifactEligible(workflow_data, kind)
}

#' Resolve coordinator-owned DA contrasts
#' @param workflow_data Mutable proteomics workflow state.
#' @param da_data Optional mutable DA module state.
#' @return A normalized contrasts table.
#' @noRd
protDaResolveContrasts <- function(workflow_data, da_data = NULL) {
    if (!protDaArtifactWorkflow(workflow_data)) {
        return(protDiaDaResolveContrasts(workflow_data, da_data))
    }
    contrasts <- resolveProtContextDependency(
        "contrasts",
        workflow_data = workflow_data,
        required = TRUE
    )$value
    if (is.null(contrasts) || !is.data.frame(contrasts) ||
        nrow(contrasts) == 0L) {
        protDiaDaArtifactAbort(
            "artifact proteomics DA has no coordinator-owned contrasts",
            "multischolar_missing_prot_da_contrasts"
        )
    }
    normaliseProtDaContrastsTable(contrasts)
}

#' Resolve coordinator-owned DA annotations
#' @param workflow_data Mutable proteomics workflow state.
#' @return The current annotations table, or `NULL`.
#' @noRd
protDaResolveAnnotations <- function(workflow_data) {
    if (!protDaArtifactWorkflow(workflow_data)) {
        return(protDiaDaResolveAnnotations(workflow_data))
    }
    resolveProtContextDependency(
        "annotations",
        workflow_data = workflow_data
    )$value
}

#' Build normalized DA parameters
#' @param ... Arguments forwarded to [protDiaDaParameters()].
#' @return Safe normalized DA parameters.
#' @noRd
protDaParameters <- function(...) protDiaDaParameters(...)

#' Prepare a descriptor-dispatched DA artifact run
#' @param workflow_data Mutable proteomics workflow state.
#' @param contrasts Normalized contrasts table.
#' @param results Current in-memory DA result list.
#' @param parameters Normalized scientific parameters.
#' @param failure_injector Optional failure injector used by tests.
#' @param now Run timestamp.
#' @return A prepared immutable DA run.
#' @noRd
prepareProtDaArtifactRun <- function(
    workflow_data,
    contrasts,
    results,
    parameters,
    failure_injector = NULL,
    now = Sys.time()
) {
    prepareProtDiaDaArtifactRun(
        workflow_data,
        contrasts,
        results,
        parameters,
        failure_injector,
        now,
        eligibility_fn = protDaArtifactEligible
    )
}

#' Publish a descriptor-dispatched DA artifact run
#' @param workflow_data Mutable proteomics workflow state.
#' @param prepared Prepared DA run.
#' @param failure_injector Optional failure injector used by tests.
#' @return Published DA run and current index.
#' @noRd
publishProtDaArtifactRun <- function(
    workflow_data,
    prepared,
    failure_injector = NULL
) {
    publishProtDiaDaArtifactRun(workflow_data, prepared, failure_injector)
}

#' Restore a descriptor-dispatched current DA artifact index
#' @param workflow_data Mutable proteomics workflow state.
#' @return A current DA index, or `NULL`.
#' @noRd
restoreProtDaArtifactIndex <- function(workflow_data) {
    restoreProtDiaDaArtifactIndex(
        workflow_data,
        eligibility_fn = protDaArtifactEligible
    )
}

#' Query one descriptor-dispatched bounded DA page
#' @param ... Arguments forwarded to [queryProtDiaDaPage()].
#' @return A bounded keyset page.
#' @noRd
queryProtDaPage <- function(...) {
    arguments <- list(...)
    workflow_data <- arguments$workflow_data
    if (is.null(workflow_data) && length(arguments) > 0L) {
        workflow_data <- arguments[[1L]]
    }
    if (protNonDiaDaArtifactWorkflow(workflow_data) &&
        !protNonDiaDaArtifactEligible(workflow_data, "queries")) {
        protDiaDaArtifactAbort(
            "non-DIA DA server-side queries are disabled",
            "multischolar_prot_nondia_da_queries_disabled"
        )
    }
    queryProtDiaDaPage(..., eligibility_fn = protDaArtifactEligible)
}

#' Collect one bounded complete selected DA table
#' @param ... Arguments forwarded to [protDiaDaCompleteSelectedTable()].
#' @return A complete bounded selected contrast table.
#' @noRd
protDaCompleteSelectedTable <- function(...) {
    protDiaDaCompleteSelectedTable(
        ...,
        eligibility_fn = protDaArtifactEligible
    )
}

#' Resolve volcano or heatmap results by exact descriptor
#' @param ... Arguments forwarded to [protDiaDaSelectedResults()].
#' @return Current selected results and S4 object.
#' @noRd
protDaSelectedResults <- function(...) protDiaDaSelectedResults(...)

#' Resolve bounded DA table filters
#' @param ... Arguments forwarded to [protDiaDaTableFilters()].
#' @return A bounded filter list.
#' @noRd
protDaTableFilters <- function(...) protDiaDaTableFilters(...)

#' Summarize one descriptor-dispatched artifact contrast
#' @param ... Arguments forwarded to [protDiaDaArtifactSummary()].
#' @return Contrast summary counts and thresholds.
#' @noRd
protDaArtifactSummary <- function(...) protDiaDaArtifactSummary(...)

#' Test any descriptor-dispatched proteomics DA artifact index
#' @param value Candidate DA result index.
#' @return A scalar logical.
#' @noRd
isProtDaArtifactIndex <- function(value) isProtDiaDaArtifactIndex(value)

#' Build exact session stage evidence for DA source binding
#' @param workflow_data Mutable proteomics workflow state.
#' @return Import and design stage evidence.
#' @noRd
protDaSessionStageEvidence <- function(workflow_data) {
    if (protNonDiaDaArtifactWorkflow(workflow_data)) {
        descriptor <- protNonDiaNormDescriptor(workflow_data)
        return(protNonDiaSessionStageEvidence(workflow_data, descriptor))
    }
    protDiaSessionStageEvidence(workflow_data)
}

#' Validate import/design DA source refs by exact workflow slug
#' @param store Validated artifact store.
#' @param stage Session stage evidence.
#' @param roles Exact expected roles.
#' @param stage_id Stage identifier.
#' @return Validated stage evidence.
#' @noRd
protDaSessionValidateStageForStore <- function(
    store,
    stage,
    roles,
    stage_id
) {
    if (identical(store$labels$workflow_slug, "prot_dia")) {
        return(protDiaSessionValidateStage(store, stage, roles, stage_id))
    }
    matches <- Filter(
        \(descriptor) identical(
            descriptor$identity$workflow_slug,
            store$labels$workflow_slug
        ),
        protNonDiaReadthroughDescriptors()
    )
    if (length(matches) != 1L) {
        protDiaDaArtifactAbort(
            "DA source stage has no exact workflow descriptor",
            "multischolar_invalid_prot_da_source_descriptor"
        )
    }
    protNonDiaSessionValidateStage(
        store,
        stage,
        roles,
        stage_id,
        protNonDiaArtifactReadthroughAdapter(matches[[1L]])
    )
}

#' Apply a descriptor-dispatched current DA artifact index
#' @param workflow_data Mutable proteomics workflow state.
#' @param da_data Mutable DA module state.
#' @return A scalar logical indicating whether an index was applied.
#' @noRd
applyProtDaArtifactIndex <- function(workflow_data, da_data) {
    index <- restoreProtDaArtifactIndex(workflow_data)
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
