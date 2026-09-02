# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

#' Resolve an injected metabolomics pre-ingress receipt envelope
#'
#' @return An injected envelope or `NULL` for installed legacy policy.
#' @noRd
metabPreIngressInjectedEnvelope <- function() {
    getOption("MultiScholaR.metabolomics.preingress_envelope")
}

#' Build one exact custom metabolomics pre-ingress contract
#'
#' @param context Unbound metabolomics workflow context.
#' @param headers Bounded source headers.
#'
#' @return A validated pre-ingress contract.
#' @noRd
metabPreIngressContract <- function(context, headers) {
    if (!inherits(context, "WorkflowContext") || context$isBound() ||
        !is.character(headers) || !length(headers) || anyDuplicated(headers)) {
        workflowPreIngressAbort(
            "metabolomics pre-ingress context or headers are invalid",
            "multischolar_invalid_metabolomics_preingress"
        )
    }
    descriptor <- artifactMetabolomicsWorkflowDescriptor()
    capability_id <- descriptor$descriptor_id
    static_identity <- context$getStaticIdentity()
    newWorkflowPreIngressContract(
        contract_id = "metabolomics.custom.tabular.preingress.v1",
        capability_ids = capability_id,
        member_resolver = \(input) {
            lapply(seq_along(input$paths), \(index) {
                list(
                    member_id = paste0("assay_", sprintf("%04d", index)),
                    path = input$paths[[index]],
                    container_type = "plain",
                    semantic_order = as.integer(index)
                )
            })
        },
        identity_probe = \(members) {
            first <- members[[1L]]
            source <- first$canonical_path
            connection <- file(source, open = "rb")
            on.exit(close(connection), add = TRUE)
            header <- readLines(connection, n = 1L, warn = FALSE)
            extension <- tolower(tools::file_ext(source))
            identity <- c(
                static_identity,
                descriptor$identity[setdiff(
                    names(descriptor$identity),
                    names(static_identity)
                )]
            )
            list(
                identity = identity,
                capability_id = capability_id,
                capability_version = descriptor$descriptor_version,
                descriptor_digest = descriptor$descriptor_digest,
                bytes_read = if (length(header)) {
                    nchar(header[[1L]], type = "bytes")
                } else {
                    0
                },
                lines_read = as.integer(length(header)),
                schema_valid = extension %in% c("csv", "tsv", "txt") &&
                    length(headers) > 1L && all(nzchar(headers)),
                ambiguous = FALSE,
                complete_payload_materialized = FALSE
            )
        },
        max_probe_bytes = 1024^2,
        max_probe_lines = 1L
    )
}

#' Resolve custom metabolomics routing before the full reader
#'
#' @param context Unbound workflow context.
#' @param source_path Existing custom source.
#' @param format Exact resolved format.
#' @param headers Bounded source headers.
#' @param injected_envelope Optional exact receipt envelope.
#'
#' @return Routing evidence and optional bound outcome.
#' @noRd
resolveMetabPreIngress <- function(
    context,
    source_path,
    format,
    headers,
    injected_envelope = metabPreIngressInjectedEnvelope()
) {
    if (!inherits(context, "WorkflowContext") || context$isBound() ||
        !identical(format, "custom")) {
        return(list(
            status = "not_applicable",
            defer_full_import = FALSE,
            outcome = NULL
        ))
    }
    policy <- context$getStoragePolicy()
    if (identical(policy$requested_backend, "artifact")) {
        return(list(
            status = "explicit_artifact",
            defer_full_import = TRUE,
            outcome = NULL
        ))
    }
    contract <- metabPreIngressContract(context, headers)
    if (identical(policy$requested_backend, "memory")) {
        outcome <- workflowPreIngressResolve(
            contract,
            list(paths = source_path),
            dirname(normalizePath(source_path)),
            requested_backend = "memory"
        )
        return(list(
            status = outcome$status,
            defer_full_import = FALSE,
            outcome = outcome
        ))
    }
    if (is.null(injected_envelope)) {
        return(list(
            status = "installed_legacy_policy_deferred",
            defer_full_import = FALSE,
            outcome = NULL
        ))
    }
    outcome <- workflowPreIngressResolve(
        contract,
        list(paths = source_path),
        dirname(normalizePath(source_path)),
        requested_backend = "auto",
        injected_envelope = injected_envelope
    )
    if (identical(outcome$status, "bound")) {
        workflowPreIngressBindContext(context, outcome, buildArtifactPaths)
    }
    list(
        status = outcome$status,
        defer_full_import = identical(outcome$effective_backend, "artifact"),
        outcome = outcome
    )
}
