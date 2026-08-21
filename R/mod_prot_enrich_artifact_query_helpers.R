# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

protDiaEnrichQuerySpecification <- function(manifest, table) {
    filters <- list()
    if ("directionality" %in% names(table)) {
        filters$direction <- list(
            column = "directionality",
            type = "character",
            operators = "in"
        )
    }
    newArtifactQueryPageSpecification(
        query_id = paste0("proteomics.diann.enrichment.", manifest$run_id, ".v1"),
        state_role = manifest$result_table$ref$logical_key$state_role,
        projections = c(names(table), .PROT_DIA_ENRICH_ROW_ORDER_COLUMN),
        filters = filters,
        sorts = list(row_order = list(
            column = .PROT_DIA_ENRICH_ROW_ORDER_COLUMN,
            transform = "identity",
            directions = "asc"
        )),
        default_sort = list(sort_id = "row_order", direction = "asc"),
        max_rows = .PROT_DIA_ENRICH_MAX_RESULT_ROWS,
        max_bytes = .PROT_DIA_ENRICH_MAX_RESULT_BYTES
    )
}

protDiaEnrichRunIndex <- function(manifest, context, interactive_plots = NULL) {
    store <- protDiaEnrichStore(context)
    table <- protDiaEnrichReadTable(store, manifest$result_table)
    structure(
        list(
            schema = .PROT_DIA_ENRICH_INDEX_SCHEMA,
            schema_version = .PROT_DIA_ENRICH_INDEX_VERSION,
            backend = "artifact",
            run_id = manifest$run_id,
            source = manifest$source,
            parameters = manifest$parameters,
            manifest_relative_path = protDiaEnrichPaths(
                context,
                run_id = manifest$run_id
            )$run_manifest,
            manifest_digest = manifest$manifest_digest,
            result_table = manifest$result_table,
            query_specification = protDiaEnrichQuerySpecification(
                manifest,
                table
            ),
            products = manifest$products,
            requests = manifest$requests,
            interactive_plots = interactive_plots
        ),
        class = c("MultiScholaRProtDiaEnrichIndex", "list")
    )
}

protDiaEnrichPlotProductName <- function(index, direction) {
    contrast <- gsub(
        "[^A-Za-z0-9_.-]",
        "_",
        index$source$contrast_name
    )
    paste0(contrast, "_", direction, "_enrichment_results.tsv")
}

protDiaEnrichPlotProduct <- function(index, direction) {
    expected <- protDiaEnrichPlotProductName(index, direction)
    matches <- Filter(
        \(product) identical(product$name, expected),
        index$products
    )
    if (length(matches) == 0L) return(NULL)
    if (length(matches) != 1L) {
        protDiaEnrichArtifactAbort(
            "enrichment plot data product is ambiguous",
            "multischolar_prot_dia_enrichment_product_mismatch"
        )
    }
    matches[[1L]]
}

protDiaEnrichPlotRecord <- function(index, direction) {
    matches <- Filter(
        \(record) identical(record$direction, direction),
        index$requests
    )
    if (length(matches) != 1L) {
        protDiaEnrichArtifactAbort(
            "enrichment plot request provenance is incomplete",
            "multischolar_invalid_prot_dia_enrichment_manifest"
        )
    }
    matches[[1L]]
}

protDiaEnrichReadPlotData <- function(workflow_data, index, direction) {
    record <- protDiaEnrichPlotRecord(index, direction)
    if (!identical(record$status, "succeeded") || record$response$rows == 0L) {
        return(NULL)
    }
    product <- protDiaEnrichPlotProduct(index, direction)
    if (is.null(product) || product$bytes > .PROT_DIA_ENRICH_MAX_RESULT_BYTES) {
        protDiaEnrichArtifactAbort(
            "enrichment plot data exceeds its bounded product contract",
            "multischolar_prot_dia_enrichment_complete_view_limit"
        )
    }
    context <- workflow_data$workflow_context
    path <- artifactResolveContainedPath(
        context$getProjectRoot(),
        product$relative_path,
        must_exist = TRUE
    )
    plot_data <- tryCatch(
        readr::read_tsv(
            path,
            show_col_types = FALSE,
            progress = FALSE,
            name_repair = "minimal"
        ) |>
            as.data.frame(stringsAsFactors = FALSE),
        error = \(error) protDiaEnrichArtifactAbort(
            "enrichment plot data could not be read",
            "multischolar_corrupt_prot_dia_enrichment_product",
            parent = error
        )
    )
    if (nrow(plot_data) != record$response$rows ||
        nrow(plot_data) > .PROT_DIA_ENRICH_MAX_RESULT_ROWS) {
        protDiaEnrichArtifactAbort(
            "enrichment plot data differs from its response provenance",
            "multischolar_prot_dia_enrichment_product_mismatch"
        )
    }
    plot_data
}

protDiaEnrichBuildRestoredPlot <- function(index, plot_data, direction) {
    if (is.null(plot_data)) return(NULL)
    plot <- if (identical(index$parameters$backend, "gprofiler2")) {
        result <- list(
            result = plot_data,
            meta = list(query_metadata = list(
                user_threshold = index$parameters$q_cutoff
            ))
        )
        buildGprofilerEnrichmentPlots(
            result,
            index$source$contrast_name,
            direction
        )$interactive
    } else {
        buildClusterProfilerEnrichmentPlots(
            plot_data,
            index$source$contrast_name,
            direction,
            index$parameters$q_cutoff
        )$interactive
    }
    plot
}

protDiaEnrichRestorePlots <- function(workflow_data, index) {
    plots <- lapply(c("up", "down"), \(direction) {
        plot_data <- protDiaEnrichReadPlotData(
            workflow_data,
            index,
            direction
        )
        withCallingHandlers(
            protDiaEnrichBuildRestoredPlot(index, plot_data, direction),
            warning = \(warning) {
                unsupported <- startsWith(
                    conditionMessage(warning),
                    "geom_GeomTextRepel() has yet to be implemented"
                )
                if (isTRUE(unsupported)) invokeRestart("muffleWarning")
            }
        )
    })
    stats::setNames(plots, c("up", "down"))
}

isProtDiaEnrichArtifactIndex <- function(value) {
    inherits(value, "MultiScholaRProtDiaEnrichIndex") &&
        identical(value$schema, .PROT_DIA_ENRICH_INDEX_SCHEMA) &&
        identical(value$schema_version, .PROT_DIA_ENRICH_INDEX_VERSION) &&
        identical(value$backend, "artifact")
}

protDiaEnrichArtifactContrastChoices <- function(index) {
    if (!isProtDiaDaArtifactIndex(index)) {
        protDiaEnrichArtifactAbort(
            "enrichment contrast choices require a DIA-NN DA artifact index",
            "multischolar_invalid_prot_dia_enrichment_da_index"
        )
    }
    values <- vapply(index$contrasts, function(entry) {
        if (workflowCapabilityScalarString(entry$friendly_name)) {
            entry$friendly_name
        } else {
            entry$contrast_name
        }
    }, character(1))
    list(
        contrastsAvailable = values,
        contrastChoices = stats::setNames(values, values),
        source = "artifact_index"
    )
}

initialiseProtDiaEnrichArtifactData <- function(workflow_data, enrichment_data) {
    index <- protDiaEnrichDaIndex(workflow_data)
    current_state <- workflowStateCurrentName(workflow_data$state_manager)
    current_s4 <- workflow_data$state_manager$getState(current_state)
    choices <- protDiaEnrichArtifactContrastChoices(index)
    enrichment_data$current_s4_object <- current_s4
    enrichment_data$da_artifact_index <- index
    enrichment_data$da_results_data <- index
    enrichment_data$contrasts_available <- choices$contrastsAvailable

    enrichment_index <- restoreProtDiaEnrichArtifactIndex(workflow_data)
    if (!is.null(enrichment_index)) {
        payload <- protDiaEnrichReactivePayload(workflow_data, enrichment_index)
        enrichment_data$gprofiler_results <- payload$gprofiler_results
        enrichment_data$clusterprofiler_results <-
            payload$clusterprofiler_results
        enrichment_data$stringdb_results <- NULL
        enrichment_data$enrichment_results_full <- enrichment_index
        enrichment_data$analysis_complete <- TRUE
        workflow_data$enrichment_analysis_results <- payload
        workflow_data$enrichment_artifact_index <- enrichment_index
    }
    list(
        index = index,
        currentS4 = current_s4,
        contrastConfig = choices,
        enrichmentIndex = enrichment_index
    )
}

restoreProtDiaEnrichArtifactIndex <- function(workflow_data) {
    if (!protDiaEnrichArtifactEligible(workflow_data, "queries")) return(NULL)
    context <- workflow_data$workflow_context
    current <- artifactResolveContainedPath(
        context$getProjectRoot(), protDiaEnrichPaths(context)$current
    )
    if (!file.exists(current)) return(NULL)
    manifest <- protDiaEnrichReadJson(
        current,
        function(value) protDiaEnrichValidateRunManifest(value, context)
    )
    da_index <- protDiaEnrichDaIndex(workflow_data)
    entry <- protDiaDaIndexEntry(
        da_index,
        manifest$source$contrast_name
    )
    valid <- identical(manifest$source$da_run_id, da_index$run_id) &&
        identical(manifest$source$da_manifest_digest, da_index$manifest_digest) &&
        identical(manifest$source$contrast_id, entry$contrast_id) &&
        identical(
            manifest$source$contrast_manifest_digest,
            entry$manifest_digest
        )
    if (!isTRUE(valid)) return(NULL)
    index <- protDiaEnrichRunIndex(manifest, context)
    index$interactive_plots <- protDiaEnrichRestorePlots(workflow_data, index)
    index
}

queryProtDiaEnrichPage <- function(
    workflow_data,
    index,
    projections = NULL,
    direction = NULL,
    cursor = NULL,
    limit = 100L,
    resource_policy = NULL
) {
    if (!protDiaEnrichArtifactEligible(workflow_data, "queries") ||
        !isProtDiaEnrichArtifactIndex(index)) {
        protDiaEnrichArtifactAbort(
            "DIA-NN enrichment query is disabled or invalid",
            "multischolar_prot_dia_enrichment_queries_disabled"
        )
    }
    filters <- list()
    if (!is.null(direction)) {
        filters$direction <- list(
            operator = "in",
            value = as.character(direction)
        )
    }
    if (is.null(projections)) {
        projections <- setdiff(
            index$query_specification$projections,
            .PROT_DIA_ENRICH_ROW_ORDER_COLUMN
        )
    }
    queryArtifactPage(
        store = protDiaEnrichStore(workflow_data$workflow_context),
        ref = index$result_table$ref,
        specification = index$query_specification,
        projections = projections,
        filters = filters,
        cursor = cursor,
        limit = limit,
        resource_policy = resource_policy
    )
}

protDiaEnrichCompleteTable <- function(
    workflow_data,
    index,
    resource_policy = NULL
) {
    rows <- as.integer(index$result_table$rows)
    page <- queryProtDiaEnrichPage(
        workflow_data,
        index,
        cursor = NULL,
        limit = max(1L, rows),
        resource_policy = resource_policy
    )
    if (isTRUE(page$has_more)) {
        protDiaEnrichArtifactAbort(
            "enrichment results exceed the complete interactive bound",
            "multischolar_prot_dia_enrichment_complete_view_limit"
        )
    }
    protDiaEnrichRestoreTable(
        page$data,
        index$result_table$list_columns,
        index$result_table$list_column_classes
    )
}

protDiaEnrichIndexWithPlots <- function(index, enrichment_results) {
    if (!isProtDiaEnrichArtifactIndex(index)) return(index)
    contrast <- index$source$contrast_name
    plots <- tryCatch(
        enrichment_results@enrichment_plotly[[contrast]],
        error = function(error) NULL
    )
    index$interactive_plots <- plots
    index
}

protDiaEnrichInteractivePlots <- function(value, contrast) {
    if (isProtDiaEnrichArtifactIndex(value)) {
        return(value$interactive_plots)
    }
    value@enrichment_plotly[[contrast]]
}

protDiaEnrichReactivePayload <- function(workflow_data, index) {
    results <- protDiaEnrichCompleteTable(workflow_data, index)
    method <- index$parameters$backend
    list(
        gprofiler_results = if (identical(method, "gprofiler2")) {
            results
        } else {
            NULL
        },
        clusterprofiler_results = if (identical(method, "clusterprofiler")) {
            results
        } else {
            NULL
        },
        stringdb_results = NULL,
        full_enrichment_results = index,
        selected_contrast = index$source$contrast_name,
        analysis_method = method,
        parameters = list(
            up_cutoff = index$parameters$up_cutoff,
            down_cutoff = index$parameters$down_cutoff,
            q_cutoff = index$parameters$q_cutoff,
            organism_taxid = index$parameters$organism_taxid
        )
    )
}

finalizeProtDiaEnrichArtifactResults <- function(
    selectedContrast,
    setupConfig,
    enrichmentResults,
    enrichmentData,
    workflowData,
    input,
    methodInfo,
    executionContext,
    capturePostProcessResultsFn = captureProtEnrichPostProcessResults,
    buildAllContrastResultsFn = buildProtEnrichAllContrastResults,
    resolveSelectedContrastResultsFn = resolveProtEnrichSelectedContrastResults,
    propagateUiParamsFn = propagateProtEnrichUiParams,
    updateStateManagerUiParamsFn = updateProtEnrichStateManagerUiParams,
    completeTabStatusFn = completeProtEnrichTabStatus,
    completeProgressFn = completeProtEnrichProgress,
    failure_injector = NULL,
    catFn = cat
) {
    capturePostProcessResultsFn(
        selectedContrast = selectedContrast,
        enrichmentResults = enrichmentResults,
        enrichmentData = enrichmentData,
        contrastsTbl = setupConfig$contrastsTbl,
        methodInfo = methodInfo,
        buildAllContrastResultsFn = buildAllContrastResultsFn,
        resolveSelectedContrastResultsFn = resolveSelectedContrastResultsFn,
        catFn = catFn
    )
    records <- protDiaEnrichValidateExecutionRecords(
        executionContext,
        methodInfo$method,
        setupConfig$rawContrastName
    )
    result_table <- if (identical(methodInfo$method, "gprofiler2")) {
        enrichmentData$gprofiler_results
    } else {
        enrichmentData$clusterprofiler_results
    }
    if (is.null(result_table) || ncol(result_table) == 0L) {
        result_table <- data.frame(
            directionality = character(),
            analysis_method = character(),
            stringsAsFactors = FALSE
        )
    }
    parameters <- protDiaEnrichParameters(
        input,
        methodInfo,
        filter_applied = setupConfig$organismFilterApplied,
        filter_stats = setupConfig$filterStats
    )
    prepared <- protDiaEnrichPersistRun(
        workflowData,
        setupConfig$artifactSelection,
        parameters,
        records,
        result_table,
        setupConfig$pathwayDir,
        failure_injector
    )
    prepared <- publishProtDiaEnrichRun(
        workflowData,
        prepared,
        failure_injector
    )
    index <- protDiaEnrichRunIndex(
        prepared$manifest,
        workflowData$workflow_context
    )
    index <- protDiaEnrichIndexWithPlots(index, enrichmentResults)
    payload <- protDiaEnrichReactivePayload(workflowData, index)
    enrichmentData$gprofiler_results <- payload$gprofiler_results
    enrichmentData$clusterprofiler_results <- payload$clusterprofiler_results
    enrichmentData$stringdb_results <- NULL
    enrichmentData$all_enrichment_results <- list(
        selected = payload[c(
            "gprofiler_results",
            "clusterprofiler_results",
            "stringdb_results"
        )]
    )
    enrichmentData$enrichment_results_full <- index
    enrichmentData$analysis_complete <- TRUE
    workflowData$enrichment_analysis_results <- payload
    workflowData$enrichment_artifact_index <- index

    ui <- propagateUiParamsFn(
        currentS4Object = enrichmentData$current_s4_object,
        workflowData = workflowData,
        selectedContrast = selectedContrast,
        methodInfo = methodInfo,
        upCutoff = input$up_cutoff,
        downCutoff = input$down_cutoff,
        qCutoff = input$q_cutoff,
        organismTaxid = input$organism_taxid
    )
    enrichmentData$current_s4_object <- ui$currentS4Object
    updateStateManagerUiParamsFn(
        workflowData = workflowData,
        storedUiParams = ui$storedUiParams
    )
    completeTabStatusFn(workflowData = workflowData)
    completeProgressFn()
    list(
        selectedContrast = selectedContrast,
        rawContrastName = setupConfig$rawContrastName,
        analysisMethod = methodInfo$method,
        analysisComplete = TRUE,
        organismFilterApplied = setupConfig$organismFilterApplied,
        filterStats = setupConfig$filterStats,
        enrichmentResults = index
    )
}
