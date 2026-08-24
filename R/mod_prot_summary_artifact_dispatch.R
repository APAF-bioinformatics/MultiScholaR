# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

#' Resolve an exact non-DIA summary artifact switch
#'
#' @param workflow_data Mutable proteomics workflow state.
#' @param kind Report dependency reads or final S4 export.
#'
#' @return One of `"enabled"` or `"disabled"`.
#' @noRd
protNonDiaSummaryMode <- function(
    workflow_data,
    kind = c("report_reads", "final_export")
) {
    kind <- match.arg(kind)
    suffix <- switch(
        kind,
        report_reads = "summary_artifact_reads",
        final_export = "summary_final_export"
    )
    value <- getOption(paste0("multischolar.prot_nondia.", suffix), "enabled")
    descriptor <- protNonDiaNormDescriptor(workflow_data)
    if (!is.null(descriptor)) {
        option <- paste0(
            "multischolar.",
            descriptor$identity$workflow_slug,
            ".",
            suffix
        )
        value <- getOption(option, value)
    }
    match.arg(value, c("enabled", "disabled"))
}

#' Test exact non-DIA summary artifact ownership
#'
#' @param workflow_data Mutable proteomics workflow state.
#'
#' @return A scalar logical.
#' @noRd
protNonDiaSummaryArtifactWorkflow <- function(workflow_data) {
    protNonDiaDaArtifactWorkflow(workflow_data)
}

#' Test exact non-DIA summary artifact eligibility
#'
#' @param workflow_data Mutable proteomics workflow state.
#' @param kind Report dependency reads or final S4 export.
#'
#' @return A scalar logical.
#' @noRd
protNonDiaSummaryArtifactEligible <- function(
    workflow_data,
    kind = c("report_reads", "final_export")
) {
    kind <- match.arg(kind)
    identical(protNonDiaSummaryMode(workflow_data, kind), "enabled") &&
        protNonDiaSummaryArtifactWorkflow(workflow_data)
}

#' Test descriptor-dispatched summary artifact eligibility
#'
#' @param workflow_data Mutable proteomics workflow state.
#' @param kind Report dependency reads or final S4 export.
#'
#' @return A scalar logical.
#' @noRd
protSummaryArtifactEligible <- function(
    workflow_data,
    kind = c("report_reads", "final_export")
) {
    kind <- match.arg(kind)
    if (protDiaDaArtifactWorkflow(workflow_data)) {
        return(protDiaSummaryArtifactEligible(workflow_data, kind))
    }
    protNonDiaSummaryArtifactEligible(workflow_data, kind)
}

#' Resolve artifact config without hydrating the scientific S4 object
#'
#' @param workflow_data Mutable proteomics workflow state.
#'
#' @return Final state name and generation-owned config, or `NULL`.
#' @noRd
protSummaryArtifactConfig <- function(workflow_data) {
    if (!protSummaryArtifactEligible(workflow_data, "report_reads")) {
        return(NULL)
    }
    manager <- workflow_data$state_manager
    state_name <- protDiaSummaryFinalStateName(manager)
    metadata <- manager$getStateMetadata(state_name)
    if (!is.list(metadata$config)) {
        protDiaSummaryAbort(
            "summary final generation has no explicit config",
            "multischolar_missing_prot_dia_summary_dependency",
            dependency = "config_list"
        )
    }
    list(
        stateName = state_name,
        configList = metadata$config,
        workflowType = workflow_data$workflow_context$getIdentity()$workflow_type
    )
}

#' Build a descriptor-dispatched summary session audit
#'
#' @param workflow_data Mutable proteomics workflow state.
#'
#' @return A bounded summary audit, or `NULL`.
#' @noRd
protSummarySessionAudit <- function(workflow_data) {
    protDiaSummarySessionAudit(
        workflow_data,
        eligibility_fn = protSummaryArtifactEligible
    )
}

#' Prepare descriptor-dispatched summary dependencies
#'
#' @param ... Arguments forwarded to [prepareProtDiaSummaryDependencies()].
#'
#' @return Bounded report dependencies or an exact final S4 dependency bundle.
#' @noRd
prepareProtSummaryDependencies <- function(...) {
    prepareProtDiaSummaryDependencies(
        ...,
        eligibility_fn = protSummaryArtifactEligible
    )
}

#' Release descriptor-dispatched summary dependencies
#'
#' @param dependencies Summary dependency bundle.
#'
#' @return A scalar logical, invisibly.
#' @noRd
releaseProtSummaryDependencies <- function(dependencies) {
    releaseProtDiaSummaryDependencies(dependencies)
}

#' Prepare declared report dependencies
#'
#' @param dependencies Summary dependency bundle.
#' @param template_path Project-contained report template path.
#'
#' @return Summary dependencies with a bounded report manifest.
#' @noRd
prepareProtReportDependencies <- function(dependencies, template_path) {
    prepareProtDiaReportDependencies(dependencies, template_path)
}

#' Write an exact descriptor-dispatched final S4 export
#'
#' @param ... Arguments forwarded to [writeProtDiaSummaryFinalExport()].
#'
#' @return Final export product and metadata paths.
#' @noRd
writeProtSummaryFinalExport <- function(...) {
    writeProtDiaSummaryFinalExport(...)
}

#' Record a descriptor-dispatched report product
#'
#' @param dependencies Summary dependency bundle.
#' @param path Rendered report path.
#'
#' @return Report product and metadata paths.
#' @noRd
recordProtSummaryReportProduct <- function(dependencies, path) {
    recordProtDiaSummaryReportProduct(dependencies, path)
}
