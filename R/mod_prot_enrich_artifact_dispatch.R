# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

#' Resolve an exact non-DIA enrichment artifact switch
#'
#' @param workflow_data Mutable proteomics workflow state.
#' @param kind Persistence or server-side query direction.
#'
#' @return One of `"enabled"` or `"disabled"`.
#' @noRd
protNonDiaEnrichArtifactMode <- function(
    workflow_data,
    kind = c("persistence", "queries")
) {
    kind <- match.arg(kind)
    descriptor <- protNonDiaNormDescriptor(workflow_data)
    generic_option <- paste0(
        "multischolar.prot_nondia.enrichment_",
        kind
    )
    value <- getOption(generic_option, "enabled")
    if (!is.null(descriptor)) {
        tuple_option <- paste0(
            "multischolar.",
            descriptor$identity$workflow_slug,
            ".enrichment_",
            kind
        )
        value <- getOption(tuple_option, value)
    }
    match.arg(value, c("enabled", "disabled"))
}

#' Test exact non-DIA enrichment artifact ownership
#'
#' @param workflow_data Mutable proteomics workflow state.
#'
#' @return A scalar logical.
#' @noRd
protNonDiaEnrichArtifactWorkflow <- function(workflow_data) {
    protNonDiaDaArtifactWorkflow(workflow_data)
}

#' Test exact non-DIA enrichment artifact eligibility
#'
#' @param workflow_data Mutable proteomics workflow state.
#' @param kind Persistence or server-side query direction.
#'
#' @return A scalar logical.
#' @noRd
protNonDiaEnrichArtifactEligible <- function(
    workflow_data,
    kind = "persistence"
) {
    identical(
        protNonDiaEnrichArtifactMode(workflow_data, kind),
        "enabled"
    ) && protNonDiaEnrichArtifactWorkflow(workflow_data)
}

#' Test any supported proteomics enrichment artifact workflow
#'
#' @param workflow_data Mutable proteomics workflow state.
#'
#' @return A scalar logical.
#' @noRd
protEnrichArtifactWorkflow <- function(workflow_data) {
    protDiaDaArtifactWorkflow(workflow_data) ||
        protNonDiaEnrichArtifactWorkflow(workflow_data)
}

#' Test descriptor-dispatched enrichment artifact eligibility
#'
#' @param workflow_data Mutable proteomics workflow state.
#' @param kind Persistence or server-side query direction.
#'
#' @return A scalar logical.
#' @noRd
protEnrichArtifactEligible <- function(
    workflow_data,
    kind = "persistence"
) {
    if (protDiaDaArtifactWorkflow(workflow_data)) {
        return(protDiaEnrichArtifactEligible(workflow_data, kind))
    }
    protNonDiaEnrichArtifactEligible(workflow_data, kind)
}

#' Resolve the descriptor-dispatched DA source index
#'
#' @param workflow_data Mutable proteomics workflow state.
#'
#' @return A validated DA artifact index.
#' @noRd
protEnrichDaIndex <- function(workflow_data) {
    protDiaEnrichDaIndex(
        workflow_data,
        restore_index_fn = restoreProtDaArtifactIndex,
        index_predicate = isProtDaArtifactIndex
    )
}

#' Resolve one descriptor-dispatched DA contrast
#'
#' @param workflow_data Mutable proteomics workflow state.
#' @param selected_contrast Selected raw, friendly, or full contrast name.
#'
#' @return The current DA index and exact contrast entry.
#' @noRd
protEnrichSelectedContrast <- function(workflow_data, selected_contrast) {
    protDiaEnrichSelectedContrast(
        workflow_data,
        selected_contrast,
        da_index_fn = protEnrichDaIndex
    )
}

#' Build one descriptor-dispatched enrichment input object
#'
#' @param workflow_data Mutable proteomics workflow state.
#' @param selected Exact DA index selection.
#' @param current_s4_object Optional current protein S4 object.
#'
#' @return A one-contrast `da_results_for_enrichment` object.
#' @noRd
protEnrichOneContrastObject <- function(
    workflow_data,
    selected,
    current_s4_object = NULL
) {
    protDiaEnrichOneContrastObject(
        workflow_data,
        selected,
        current_s4_object,
        complete_table_fn = protDaCompleteSelectedTable,
        resolve_contrasts_fn = protDaResolveContrasts
    )
}

#' Prepare exact descriptor-owned enrichment inputs
#'
#' @param workflow_data Mutable proteomics workflow state.
#' @param selected_contrast Selected raw, friendly, or full contrast name.
#' @param current_s4_object Optional current protein S4 object.
#'
#' @return Exact DA, design, annotation, and contrast bindings.
#' @noRd
protEnrichExplicitSetup <- function(
    workflow_data,
    selected_contrast,
    current_s4_object = NULL
) {
    protDiaEnrichExplicitSetup(
        workflow_data,
        selected_contrast,
        current_s4_object,
        selected_contrast_fn = protEnrichSelectedContrast,
        one_contrast_fn = protEnrichOneContrastObject,
        resolve_annotations_fn = protDaResolveAnnotations
    )
}

#' Prepare a descriptor-dispatched enrichment artifact analysis
#'
#' @param ... Arguments forwarded to
#'   [prepareProtDiaEnrichArtifactAnalysisSetup()].
#'
#' @return Artifact-bound enrichment analysis setup.
#' @noRd
prepareProtEnrichArtifactAnalysisSetup <- function(...) {
    prepareProtDiaEnrichArtifactAnalysisSetup(
        ...,
        explicitSetupFn = protEnrichExplicitSetup
    )
}

#' Resolve the service mode for an exact proteomics workflow
#'
#' @param workflow_data Mutable proteomics workflow state.
#'
#' @return One of `"auto"`, `"live"`, or `"replay"`.
#' @noRd
protEnrichServiceMode <- function(workflow_data) {
    if (protDiaDaArtifactWorkflow(workflow_data)) {
        return(getOption(
            "multischolar.prot_dia.enrichment_service_mode",
            "auto"
        ))
    }
    descriptor <- protNonDiaNormDescriptor(workflow_data)
    value <- getOption(
        "multischolar.prot_nondia.enrichment_service_mode",
        "auto"
    )
    if (!is.null(descriptor)) {
        option <- paste0(
            "multischolar.",
            descriptor$identity$workflow_slug,
            ".enrichment_service_mode"
        )
        value <- getOption(option, value)
    }
    match.arg(value, c("auto", "live", "replay"))
}

#' Construct a descriptor-dispatched enrichment execution context
#'
#' @param workflow_data Mutable proteomics workflow state.
#' @param mode Optional service mode override.
#' @param sleep_fn Retry wait function.
#' @param failure_injector Optional test-only failure callback.
#'
#' @return An enrichment execution provenance context.
#' @noRd
newProtEnrichExecutionContext <- function(
    workflow_data,
    mode = NULL,
    sleep_fn = Sys.sleep,
    failure_injector = NULL
) {
    if (is.null(mode)) mode <- protEnrichServiceMode(workflow_data)
    newProtDiaEnrichExecutionContext(
        workflow_data,
        mode,
        sleep_fn,
        failure_injector
    )
}

#' Build a workflow-specific enrichment artifact index
#'
#' @param manifest Validated enrichment run manifest.
#' @param context Bound workflow context.
#' @param interactive_plots Optional bounded selected plots.
#'
#' @return A bounded enrichment artifact index.
#' @noRd
protEnrichRunIndex <- function(
    manifest,
    context,
    interactive_plots = NULL
) {
    workflow_slug <- context$getIdentity()$workflow_slug
    namespace <- if (identical(workflow_slug, "prot_dia")) {
        "proteomics.diann.enrichment"
    } else {
        paste0("proteomics.", workflow_slug, ".enrichment")
    }
    protDiaEnrichRunIndex(
        manifest,
        context,
        interactive_plots,
        query_namespace = namespace
    )
}

#' Restore the current descriptor-dispatched enrichment index
#'
#' @param workflow_data Mutable proteomics workflow state.
#'
#' @return A current enrichment index, or `NULL`.
#' @noRd
restoreProtEnrichArtifactIndex <- function(workflow_data) {
    restoreProtDiaEnrichArtifactIndex(
        workflow_data,
        eligibility_fn = protEnrichArtifactEligible,
        da_index_fn = protEnrichDaIndex,
        run_index_fn = protEnrichRunIndex
    )
}

#' Query one descriptor-dispatched bounded enrichment page
#'
#' @param ... Arguments forwarded to [queryProtDiaEnrichPage()].
#'
#' @return A bounded keyset page.
#' @noRd
queryProtEnrichPage <- function(...) {
    queryProtDiaEnrichPage(
        ...,
        eligibility_fn = protEnrichArtifactEligible
    )
}

#' Collect a bounded complete enrichment table
#'
#' @param ... Arguments forwarded to [protDiaEnrichCompleteTable()].
#'
#' @return The complete selected enrichment result table.
#' @noRd
protEnrichCompleteTable <- function(...) {
    protDiaEnrichCompleteTable(
        ...,
        query_page_fn = queryProtEnrichPage
    )
}

#' Build the bounded reactive enrichment payload
#'
#' @param workflow_data Mutable proteomics workflow state.
#' @param index Validated enrichment artifact index.
#'
#' @return Bounded selected results and parameters for renderers.
#' @noRd
protEnrichReactivePayload <- function(workflow_data, index) {
    protDiaEnrichReactivePayload(
        workflow_data,
        index,
        complete_table_fn = protEnrichCompleteTable
    )
}

#' Resolve contrast choices from any proteomics DA artifact index
#'
#' @param index Validated proteomics DA artifact index.
#'
#' @return Contrast choices for the enrichment UI.
#' @noRd
protEnrichArtifactContrastChoices <- function(index) {
    protDiaEnrichArtifactContrastChoices(
        index,
        index_predicate = isProtDaArtifactIndex
    )
}

#' Initialise descriptor-dispatched enrichment reactive state
#'
#' @param workflow_data Mutable proteomics workflow state.
#' @param enrichment_data Mutable enrichment module state.
#'
#' @return Initial DA and enrichment artifact state.
#' @noRd
initialiseProtEnrichArtifactData <- function(
    workflow_data,
    enrichment_data
) {
    initialiseProtDiaEnrichArtifactData(
        workflow_data,
        enrichment_data,
        da_index_fn = protEnrichDaIndex,
        contrast_choices_fn = protEnrichArtifactContrastChoices,
        restore_enrichment_fn = restoreProtEnrichArtifactIndex,
        reactive_payload_fn = protEnrichReactivePayload
    )
}

#' Finalize descriptor-dispatched enrichment artifacts
#'
#' @param ... Arguments forwarded to
#'   [finalizeProtDiaEnrichArtifactResults()].
#'
#' @return Final analysis state and enrichment artifact index.
#' @noRd
finalizeProtEnrichArtifactResults <- function(...) {
    finalizeProtDiaEnrichArtifactResults(
        ...,
        runIndexFn = protEnrichRunIndex,
        reactivePayloadFn = protEnrichReactivePayload
    )
}

#' Test any proteomics enrichment artifact index
#'
#' @param value Candidate enrichment result value.
#'
#' @return A scalar logical.
#' @noRd
isProtEnrichArtifactIndex <- function(value) {
    isProtDiaEnrichArtifactIndex(value)
}

#' Resolve artifact or legacy in-memory interactive plots
#'
#' @param value Artifact index or legacy enrichment S4 object.
#' @param contrast Exact selected contrast name.
#'
#' @return The selected interactive enrichment plots.
#' @noRd
protEnrichInteractivePlots <- function(value, contrast) {
    protDiaEnrichInteractivePlots(value, contrast)
}
