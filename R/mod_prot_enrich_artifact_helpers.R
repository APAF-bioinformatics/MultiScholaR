# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

.PROT_DIA_ENRICH_RUN_SCHEMA <- "multischolar.prot_dia_enrichment_run"
.PROT_DIA_ENRICH_RUN_VERSION <- 1L
.PROT_DIA_ENRICH_RESPONSE_SCHEMA <- "multischolar.prot_dia_enrichment_response"
.PROT_DIA_ENRICH_RESPONSE_VERSION <- 1L
.PROT_DIA_ENRICH_INDEX_SCHEMA <- "multischolar.prot_dia_enrichment_index"
.PROT_DIA_ENRICH_INDEX_VERSION <- 1L
.PROT_DIA_ENRICH_MAX_RESULT_ROWS <- 250000L
.PROT_DIA_ENRICH_MAX_RESULT_BYTES <- 256L * 1024L * 1024L
.PROT_DIA_ENRICH_ROW_ORDER_COLUMN <- ".multischolar_enrichment_row_order"

protDiaEnrichArtifactAbort <- function(message, class, ...) {
    rlang::abort(
        message,
        class = c(class, "multischolar_prot_dia_enrichment_artifact_error"),
        ...
    )
}

protDiaEnrichSensitiveName <- function(value) {
    normalized <- gsub("[^a-z0-9]+", "_", tolower(value))
    grepl(
        "(^|_)(authorization|credential|password|passwd|secret|token|api_key|apikey)(_|$)",
        normalized
    )
}

protDiaEnrichSecretUrl <- function(value) {
    is.character(value) && any(grepl(
        paste0(
            "^[a-z][a-z0-9+.-]*://.*(",
            "[?&](authorization|credential|password|passwd|secret|token|api_?key)=",
            "|//[^/@:]+:[^/@]+@)"
        ),
        value,
        ignore.case = TRUE,
        perl = TRUE
    ))
}

protDiaEnrichAssertSafeMetadata <- function(value, owner = "metadata") {
    artifactWorkflowStateAssertSafeMetadata(value, owner)
    inspect <- function(candidate, candidate_owner) {
        if (protDiaEnrichSecretUrl(candidate)) {
            protDiaEnrichArtifactAbort(
                sprintf("%s contains a URL carrying credentials", candidate_owner),
                "multischolar_unsafe_artifact_state_metadata"
            )
        }
        if (!is.list(candidate)) return(invisible(TRUE))
        candidate_names <- names(candidate)
        if (!is.null(candidate_names) &&
            any(protDiaEnrichSensitiveName(candidate_names))) {
            protDiaEnrichArtifactAbort(
                sprintf("%s contains a credential field", candidate_owner),
                "multischolar_unsafe_artifact_state_metadata"
            )
        }
        for (index in seq_along(candidate)) {
            label <- if (!is.null(candidate_names) &&
                nzchar(candidate_names[[index]])) {
                candidate_names[[index]]
            } else {
                as.character(index)
            }
            inspect(
                candidate[[index]],
                paste0(candidate_owner, "[[", label, "]]")
            )
        }
        invisible(TRUE)
    }
    inspect(value, owner)
    invisible(TRUE)
}

protDiaEnrichArtifactMode <- function(kind = c("persistence", "queries")) {
    kind <- match.arg(kind)
    option <- paste0("multischolar.prot_dia.enrichment_", kind)
    match.arg(getOption(option, "enabled"), c("enabled", "disabled"))
}

protDiaEnrichArtifactEligible <- function(workflow_data, kind = "persistence") {
    identical(protDiaEnrichArtifactMode(kind), "enabled") &&
        protDiaDaArtifactWorkflow(workflow_data)
}

protDiaEnrichPaths <- function(context, run_id = NULL, request_id = NULL) {
    relative_root <- artifactNormalizeRelativePath(file.path(
        context$getPaths()$relative_paths$workflow_state_root,
        "enrichment"
    ))
    paths <- list(
        root = relative_root,
        current = artifactNormalizeRelativePath(file.path(
            relative_root, "current.json"
        )),
        responses = artifactNormalizeRelativePath(file.path(
            relative_root, "responses"
        ))
    )
    if (!is.null(run_id)) {
        artifactValidatePathComponent(run_id, "enrichment_run_id")
        paths$run_root <- artifactNormalizeRelativePath(file.path(
            relative_root, "runs", run_id
        ))
        paths$run_manifest <- artifactNormalizeRelativePath(file.path(
            paths$run_root, "manifest.json"
        ))
        paths$products <- artifactNormalizeRelativePath(file.path(
            paths$run_root, "products"
        ))
    }
    if (!is.null(request_id)) {
        artifactValidatePathComponent(request_id, "enrichment_request_id")
        paths$response_root <- artifactNormalizeRelativePath(file.path(
            paths$responses, request_id
        ))
        paths$response_manifest <- artifactNormalizeRelativePath(file.path(
            paths$response_root, "manifest.json"
        ))
    }
    paths
}

protDiaEnrichStore <- function(context) {
    identity <- context$getIdentity()
    newArtifactStore(context$getPaths(), identity$project_id)
}

protDiaEnrichExistingTableRef <- function(
    store,
    logical_key,
    encoded,
    provenance_ids
) {
    matches <- Filter(
        \(path) {
            sidecar <- artifactStoreReadSidecar(
                store,
                path,
                validate_payload = FALSE
            )
            identical(sidecar$artifact_ref$logical_key, logical_key)
        },
        artifactStoreSidecarPaths(store)
    )
    if (length(matches) == 0L) return(NULL)
    if (length(matches) != 1L) {
        protDiaEnrichArtifactAbort(
            "enrichment table has duplicate immutable payloads",
            "multischolar_prot_dia_enrichment_immutable_conflict"
        )
    }
    sidecar <- artifactStoreReadSidecar(
        store,
        matches[[1L]],
        validate_payload = TRUE
    )
    ref <- artifactStoreNormalizeRef(sidecar$artifact_ref)
    expected <- encoded$metadata
    valid <- identical(
        ref$hash_policy$semantic$digest,
        expected$semantic_digest
    ) && identical(ref$shape$rows, as.integer(expected$dimensions$rows)) &&
        identical(ref$shape$columns, as.integer(expected$dimensions$columns)) &&
        identical(ref$provenance_ids, as.character(provenance_ids))
    if (!isTRUE(valid)) {
        protDiaEnrichArtifactAbort(
            "existing enrichment table differs from the retried candidate",
            "multischolar_prot_dia_enrichment_immutable_conflict"
        )
    }
    ref
}

protDiaEnrichJsonDigest <- function(value) {
    candidate <- value
    candidate$manifest_digest <- NULL
    artifactSemanticDigest(candidate)
}

protDiaEnrichWriteJson <- function(
    value,
    path,
    validator,
    replace = FALSE,
    failure_injector = NULL,
    failure_stage = "before_enrichment_json_publication"
) {
    if (!is.function(validator)) {
        protDiaEnrichArtifactAbort(
            "enrichment JSON writes require a validator",
            "multischolar_invalid_prot_dia_enrichment_manifest"
        )
    }
    validator(value)
    parent <- dirname(path)
    if (!dir.exists(parent) && !dir.create(parent, recursive = TRUE)) {
        protDiaEnrichArtifactAbort(
            "enrichment manifest directory could not be created",
            "multischolar_prot_dia_enrichment_write_failed"
        )
    }
    temporary <- file.path(
        parent,
        paste0(".", basename(path), ".", artifactOpaqueId("tmp"), ".tmp")
    )
    on.exit(if (file.exists(temporary)) unlink(temporary, force = FALSE), add = TRUE)
    jsonlite::write_json(
        value,
        temporary,
        auto_unbox = TRUE,
        null = "null",
        na = "string",
        digits = NA,
        pretty = TRUE
    )
    restored <- jsonlite::read_json(temporary, simplifyVector = FALSE)
    validator(restored)
    if (is.function(failure_injector)) {
        failure_injector(failure_stage, value)
    }
    if (file.exists(path) && !isTRUE(replace)) {
        existing <- jsonlite::read_json(path, simplifyVector = FALSE)
        validator(existing)
        if (!identical(existing, restored)) {
            protDiaEnrichArtifactAbort(
                "immutable enrichment manifest already exists with other content",
                "multischolar_prot_dia_enrichment_immutable_conflict"
            )
        }
        return(invisible(path))
    }
    if (!isTRUE(file.rename(temporary, path))) {
        protDiaEnrichArtifactAbort(
            "enrichment manifest could not be atomically published",
            "multischolar_prot_dia_enrichment_publish_failed"
        )
    }
    invisible(path)
}

protDiaEnrichReadJson <- function(path, validator) {
    if (!file.exists(path)) return(NULL)
    value <- tryCatch(
        jsonlite::read_json(path, simplifyVector = FALSE),
        error = function(error) protDiaEnrichArtifactAbort(
            sprintf("enrichment manifest could not be read: %s", conditionMessage(error)),
            "multischolar_corrupt_prot_dia_enrichment_manifest",
            parent = error
        )
    )
    validator(value)
}

protDiaEnrichSoftware <- function(backend) {
    package <- if (identical(backend, "gprofiler2")) {
        "gprofiler2"
    } else {
        "clusterProfiler"
    }
    list(
        multischolar = protDiaNormPackageVersion("MultiScholaR"),
        backend = backend,
        backend_package = package,
        backend_version = protDiaNormPackageVersion(package),
        r = as.character(getRversion())
    )
}

protDiaEnrichParameters <- function(
    input,
    method_info,
    filter_applied = FALSE,
    filter_stats = list(
        proteins_before = 0L,
        proteins_after = 0L,
        proteins_removed = 0L
    )
) {
    filter_stats <- lapply(filter_stats[c(
        "proteins_before", "proteins_after", "proteins_removed"
    )], as.integer)
    parameters <- list(
        backend = as.character(method_info$method),
        organism_taxid = as.character(input$organism_taxid),
        organism_supported = isTRUE(method_info$supported),
        organism_name = as.character(method_info$species_name),
        up_cutoff = as.double(input$up_cutoff),
        down_cutoff = as.double(input$down_cutoff),
        q_cutoff = as.double(input$q_cutoff),
        correction_method = as.character(input$correction_method),
        exclude_iea = FALSE,
        organism_filter_enabled = isTRUE(input$enable_organism_filter),
        organism_filter_applied = isTRUE(filter_applied),
        organism_filter_stats = filter_stats
    )
    protDiaEnrichAssertSafeMetadata(
        parameters,
        "DIA-NN enrichment parameters"
    )
    parameters
}

protDiaEnrichValidateParameters <- function(value) {
    required <- c(
        "backend", "organism_taxid", "organism_supported", "organism_name",
        "up_cutoff", "down_cutoff", "q_cutoff", "correction_method",
        "exclude_iea", "organism_filter_enabled", "organism_filter_applied",
        "organism_filter_stats"
    )
    logical_fields <- c(
        "organism_supported", "exclude_iea", "organism_filter_enabled",
        "organism_filter_applied"
    )
    numeric_fields <- c("up_cutoff", "down_cutoff", "q_cutoff")
    valid <- is.list(value) && identical(names(value), required) &&
        value$backend %in% c("gprofiler2", "clusterprofiler") &&
        workflowCapabilityScalarString(value$organism_taxid) &&
        workflowCapabilityScalarString(value$organism_name) &&
        workflowCapabilityScalarString(value$correction_method) &&
        all(vapply(value[logical_fields], function(item) {
            length(item) == 1L && is.logical(item)
        }, logical(1))) && all(vapply(value[numeric_fields], function(item) {
            length(item) == 1L && is.numeric(item) && !is.na(item) &&
                is.finite(item) && item >= 0
        }, logical(1))) && value$q_cutoff <= 1
    expected_backend <- if (isTRUE(value$organism_supported)) {
        "gprofiler2"
    } else {
        "clusterprofiler"
    }
    valid <- valid && identical(value$backend, expected_backend) &&
        (!isTRUE(value$organism_filter_applied) ||
            isTRUE(value$organism_filter_enabled))
    stats <- value$organism_filter_stats
    stats_required <- c(
        "proteins_before", "proteins_after", "proteins_removed"
    )
    stats_valid <- is.list(stats) && identical(names(stats), stats_required) &&
        all(vapply(stats, function(item) {
            length(item) == 1L && is.numeric(item) && !is.na(item) &&
                is.finite(item) && item >= 0 && item == floor(item)
        }, logical(1)))
    if (isTRUE(stats_valid)) {
        stats_valid <- stats$proteins_after <= stats$proteins_before &&
            stats$proteins_removed ==
                stats$proteins_before - stats$proteins_after
    }
    if (!isTRUE(valid) || !isTRUE(stats_valid)) {
        protDiaEnrichArtifactAbort(
            "enrichment run parameters are malformed",
            "multischolar_invalid_prot_dia_enrichment_manifest"
        )
    }
    value[numeric_fields] <- lapply(value[numeric_fields], as.double)
    value$organism_filter_stats <- lapply(stats, as.integer)
    value
}

protDiaEnrichDaIndex <- function(
    workflow_data,
    restore_index_fn = restoreProtDiaDaArtifactIndex,
    index_predicate = isProtDiaDaArtifactIndex
) {
    index <- workflow_data$da_analysis_results_list
    if (!index_predicate(index)) {
        index <- restore_index_fn(workflow_data)
    }
    if (!index_predicate(index)) {
        protDiaEnrichArtifactAbort(
            "DIA-NN enrichment requires a current validated DA artifact index",
            "multischolar_missing_prot_dia_enrichment_da_index"
        )
    }
    generation <- workflow_data$state_manager$getCurrentGenerationId()
    if (!identical(index$source_generation_id, generation)) {
        protDiaEnrichArtifactAbort(
            "DIA-NN enrichment DA index is stale for the current generation",
            "multischolar_stale_prot_dia_enrichment_da_index"
        )
    }
    index
}

protDiaEnrichSelectedContrast <- function(
    workflow_data,
    selected_contrast,
    da_index_fn = protDiaEnrichDaIndex
) {
    index <- da_index_fn(workflow_data)
    entry <- protDiaDaIndexEntry(index, selected_contrast)
    list(index = index, entry = entry)
}

protDiaEnrichOneContrastObject <- function(
    workflow_data,
    selected,
    current_s4_object = NULL,
    complete_table_fn = protDiaDaCompleteSelectedTable,
    resolve_contrasts_fn = protDiaDaResolveContrasts
) {
    table <- complete_table_fn(
        workflow_data,
        selected$index,
        selected$entry$contrast_name
    )
    contrasts <- resolve_contrasts_fn(workflow_data)
    matches <- which(
        contrasts$contrasts == selected$entry$contrast_name |
            contrasts$full_format == selected$entry$full_format |
            contrasts$friendly_names == selected$entry$friendly_name
    )
    if (length(matches) != 1L) {
        protDiaEnrichArtifactAbort(
            "selected enrichment contrast does not bind uniquely to design",
            "multischolar_invalid_prot_dia_enrichment_contrast"
        )
    }
    contrasts <- tibble::as_tibble(contrasts[matches, , drop = FALSE])
    design_matrix <- workflow_data$design_matrix
    if (is.null(design_matrix) && !is.null(current_s4_object)) {
        design_matrix <- tryCatch(
            current_s4_object@design_matrix,
            error = function(error) NULL
        )
    }
    if (!is.data.frame(design_matrix)) {
        protDiaEnrichArtifactAbort(
            "DIA-NN enrichment requires an explicit design matrix",
            "multischolar_missing_prot_dia_enrichment_design"
        )
    }
    object <- methods::new("da_results_for_enrichment")
    object@contrasts <- contrasts
    object@design_matrix <- as.data.frame(design_matrix)
    object@da_data <- stats::setNames(
        list(table),
        selected$entry$contrast_name
    )
    object
}

protDiaEnrichExplicitSetup <- function(
    workflow_data,
    selected_contrast,
    current_s4_object = NULL,
    selected_contrast_fn = protDiaEnrichSelectedContrast,
    one_contrast_fn = protDiaEnrichOneContrastObject,
    resolve_annotations_fn = protDiaDaResolveAnnotations
) {
    selected <- selected_contrast_fn(workflow_data, selected_contrast)
    annotations <- resolve_annotations_fn(workflow_data)
    object <- one_contrast_fn(
        workflow_data,
        selected,
        current_s4_object
    )
    list(
        artifact = TRUE,
        index = selected$index,
        entry = selected$entry,
        raw_contrast_name = selected$entry$contrast_name,
        contrasts_tbl = object@contrasts,
        da_results = object,
        annotations = annotations
    )
}

prepareProtDiaEnrichArtifactAnalysisSetup <- function(
    selectedContrast,
    input,
    enrichmentData,
    workflowData,
    experimentPaths,
    resolveOutputDirectoriesFn = resolveProtEnrichOutputDirectories,
    resolveAnnotationMatchingFn = resolveProtEnrichAnnotationMatching,
    resolveOrganismMappingFn = resolveProtEnrichOrganismMapping,
    applyOrganismFilterFn = applyProtEnrichOrganismFilter,
    persistOrganismFilterMetadataFn = persistProtEnrichOrganismFilterMetadata,
    showNotificationFn = shiny::showNotification,
    explicitSetupFn = protDiaEnrichExplicitSetup,
    ...
) {
    current_s4 <- enrichmentData$current_s4_object
    if (is.null(current_s4)) {
        state <- workflowStateCurrentName(workflowData$state_manager)
        current_s4 <- workflowData$state_manager$getState(state)
        enrichmentData$current_s4_object <- current_s4
    }
    explicit <- explicitSetupFn(
        workflowData,
        selectedContrast,
        current_s4
    )
    enrichmentData$da_artifact_index <- explicit$index
    paths <- resolveOutputDirectoriesFn(experimentPaths)
    annotation_match <- resolveAnnotationMatchingFn(
        uniprotDatCln = explicit$annotations,
        daResultsForEnrichment = explicit$da_results,
        currentS4Object = current_s4
    )
    if (!is.null(annotation_match$annotationMatchResults)) {
        enrichmentData$annotation_match_results <-
            annotation_match$annotationMatchResults
    }
    object <- explicit$da_results
    filter_applied <- FALSE
    filter_stats <- list(
        proteins_before = 0,
        proteins_after = 0,
        proteins_removed = 0
    )
    if (isTRUE(input$enable_organism_filter)) {
        target <- as.character(input$organism_taxid)
        mapping <- resolveOrganismMappingFn(
            workflowData = workflowData,
            uniprotDatCln = explicit$annotations,
            targetTaxon = target
        )$organismMapping
        if (!is.null(mapping) && nrow(mapping) > 0L) {
            filtered <- applyOrganismFilterFn(
                daResultsForEnrichment = object,
                organismMapping = mapping,
                targetTaxon = target,
                currentS4Object = current_s4
            )
            object <- filtered$daResultsForEnrichment
            filter_applied <- filtered$filterApplied
            filter_stats <- filtered$filterStats
        } else {
            showNotificationFn(
                paste(
                    "Multi-species filtering enabled but no organism mapping",
                    "is available. Proceeding without filtering."
                ),
                type = "warning",
                duration = 8
            )
        }
    }
    persistOrganismFilterMetadataFn(
        workflowData = workflowData,
        organismFilterEnabled = input$enable_organism_filter,
        organismFilterApplied = filter_applied,
        targetTaxonId = input$organism_taxid,
        filterStats = filter_stats
    )
    list(
        artifact = TRUE,
        artifactSelection = list(index = explicit$index, entry = explicit$entry),
        rawContrastName = explicit$raw_contrast_name,
        contrastsTbl = explicit$contrasts_tbl,
        pathwayDir = paths$pathwayDir,
        daResultsForEnrichment = object,
        goAnnotations = explicit$annotations,
        organismFilterApplied = filter_applied,
        filterStats = filter_stats
    )
}

protDiaEnrichRequestId <- function(request) {
    paste0("request_", artifactSemanticDigest(request))
}

protDiaEnrichRunId <- function(
    source,
    parameters,
    records,
    software,
    result_digest
) {
    artifactRefValidateDigest(result_digest, "enrichment result digest")
    records <- records[order(match(
        vapply(records, `[[`, character(1), "direction"),
        c("up", "down")
    ))]
    request_provenance <- lapply(records, function(record) {
        list(
            direction = record$direction,
            request_digest = record$request_digest,
            status = record$status,
            execution_state = record$execution_state,
            attempts = record$attempts,
            response_digest = if (is.null(record$response)) {
                NULL
            } else {
                record$response$response_digest
            }
        )
    })
    value <- list(
        source = source,
        parameters = parameters,
        requests = request_provenance,
        software = software,
        result_digest = result_digest
    )
    paste0("enrichment_", artifactSemanticDigest(value))
}

protDiaEnrichValidateRunSource <- function(source, context) {
    required <- c(
        "da_run_id", "da_manifest_relative_path", "da_manifest_digest",
        "source_generation_id", "contrast_id", "contrast_name",
        "contrast_manifest_relative_path", "contrast_manifest_digest"
    )
    valid <- is.list(source) && identical(names(source), required) &&
        all(vapply(source, workflowCapabilityScalarString, logical(1)))
    if (!isTRUE(valid)) {
        protDiaEnrichArtifactAbort(
            "enrichment DA source binding is malformed",
            "multischolar_invalid_prot_dia_enrichment_source"
        )
    }
    artifactValidatePathComponent(source$da_run_id, "da_run_id")
    artifactValidatePathComponent(source$contrast_id, "da_contrast_id")
    artifactValidatePathComponent(
        source$source_generation_id,
        "source_generation_id"
    )
    source$da_manifest_relative_path <- artifactNormalizeRelativePath(
        source$da_manifest_relative_path
    )
    source$contrast_manifest_relative_path <- artifactNormalizeRelativePath(
        source$contrast_manifest_relative_path
    )
    artifactRefValidateDigest(source$da_manifest_digest, "DA manifest digest")
    artifactRefValidateDigest(
        source$contrast_manifest_digest,
        "DA contrast manifest digest"
    )
    path <- artifactResolveContainedPath(
        context$getProjectRoot(),
        source$da_manifest_relative_path,
        must_exist = TRUE
    )
    manifest <- protDiaDaReadJson(
        path,
        function(value) protDiaDaValidateRunManifest(value, context)
    )
    entries <- Filter(
        function(entry) identical(entry$contrast_id, source$contrast_id),
        manifest$contrasts
    )
    valid <- identical(manifest$run_id, source$da_run_id) &&
        identical(manifest$manifest_digest, source$da_manifest_digest) &&
        identical(manifest$source$generation_id, source$source_generation_id) &&
        length(entries) == 1L &&
        identical(entries[[1L]]$contrast_name, source$contrast_name) &&
        identical(
            entries[[1L]]$manifest_relative_path,
            source$contrast_manifest_relative_path
        ) && identical(
            entries[[1L]]$manifest_digest,
            source$contrast_manifest_digest
        )
    if (!isTRUE(valid)) {
        protDiaEnrichArtifactAbort(
            "enrichment DA source differs from its immutable manifests",
            "multischolar_prot_dia_enrichment_source_mismatch"
        )
    }
    source
}

protDiaEnrichProductFiles <- function(pathway_dir, contrast) {
    safe <- gsub("[^A-Za-z0-9_.-]", "_", contrast)
    pattern <- paste0(
        "^",
        safe,
        "_(up|down)_enrichment_(results\\.tsv|plot\\.(html|png))$"
    )
    list.files(pathway_dir, pattern = pattern, full.names = TRUE)
}

protDiaEnrichPersistProducts <- function(
    context,
    run_id,
    pathway_dir,
    contrast,
    failure_injector = NULL
) {
    files <- protDiaEnrichProductFiles(pathway_dir, contrast)
    if (length(files) == 0L) return(list())
    paths <- protDiaEnrichPaths(context, run_id = run_id)
    root <- artifactResolveContainedPath(
        context$getProjectRoot(),
        paths$products
    )
    if (!dir.exists(root) && !dir.create(root, recursive = TRUE)) {
        protDiaEnrichArtifactAbort(
            "enrichment product directory could not be created",
            "multischolar_prot_dia_enrichment_write_failed"
        )
    }
    lapply(files, function(source) {
        target <- file.path(root, basename(source))
        if (!file.exists(target)) {
            temporary <- paste0(
                target,
                ".",
                artifactOpaqueId("tmp"),
                ".tmp"
            )
            on.exit(
                if (file.exists(temporary)) unlink(temporary),
                add = TRUE
            )
            copied <- file.copy(source, temporary, overwrite = FALSE)
            if (isTRUE(copied) && is.function(failure_injector)) {
                failure_injector(
                    "before_enrichment_product_publication",
                    list(run_id = run_id, product_name = basename(target))
                )
            }
            if (!isTRUE(copied) || !file.rename(temporary, target)) {
                protDiaEnrichArtifactAbort(
                    "enrichment product could not be immutably published",
                    "multischolar_prot_dia_enrichment_write_failed"
                )
            }
        }
        list(
            name = basename(target),
            relative_path = artifactNormalizeRelativePath(file.path(
                paths$products,
                basename(target)
            )),
            byte_digest = artifactByteDigest(target),
            bytes = unname(as.numeric(file.info(target)$size))
        )
    })
}

protDiaEnrichExpectedProductNames <- function(source, records) {
    complete <- Filter(function(record) {
        identical(record$status, "succeeded") &&
            !is.null(record$response) && record$response$rows > 0L
    }, records)
    safe_contrast <- gsub("[^A-Za-z0-9_.-]", "_", source$contrast_name)
    expected <- unlist(lapply(complete, function(record) {
        paste0(
            safe_contrast,
            "_",
            record$direction,
            "_enrichment_",
            c("results.tsv", "plot.html", "plot.png")
        )
    }), use.names = FALSE)
    if (is.null(expected)) character() else expected
}

protDiaEnrichValidateProducts <- function(products, source, records, context, run_id) {
    if (!is.list(products)) {
        protDiaEnrichArtifactAbort(
            "enrichment product inventory is malformed",
            "multischolar_invalid_prot_dia_enrichment_manifest"
        )
    }
    expected_names <- sort(protDiaEnrichExpectedProductNames(source, records))
    observed_names <- vapply(products, function(product) {
        required <- c("name", "relative_path", "byte_digest", "bytes")
        valid <- is.list(product) && identical(names(product), required) &&
            workflowCapabilityScalarString(product$name) &&
            workflowCapabilityScalarString(product$relative_path) &&
            length(product$bytes) == 1L && is.numeric(product$bytes) &&
            !is.na(product$bytes) && is.finite(product$bytes) &&
            product$bytes > 0 && product$bytes == floor(product$bytes)
        if (!isTRUE(valid)) {
            protDiaEnrichArtifactAbort(
                "enrichment product entry is malformed",
                "multischolar_invalid_prot_dia_enrichment_manifest"
            )
        }
        artifactRefValidateDigest(product$byte_digest, "product byte digest")
        relative_path <- artifactNormalizeRelativePath(product$relative_path)
        expected_path <- artifactNormalizeRelativePath(file.path(
            protDiaEnrichPaths(context, run_id = run_id)$products,
            product$name
        ))
        path <- artifactResolveContainedPath(
            context$getProjectRoot(),
            relative_path,
            must_exist = TRUE
        )
        valid <- identical(relative_path, expected_path) &&
            identical(basename(relative_path), product$name) &&
            identical(unname(as.numeric(file.info(path)$size)),
                      as.numeric(product$bytes)) &&
            identical(artifactByteDigest(path), product$byte_digest)
        if (!isTRUE(valid)) {
            protDiaEnrichArtifactAbort(
                "enrichment product differs from its immutable manifest",
                "multischolar_prot_dia_enrichment_product_mismatch"
            )
        }
        product$name
    }, character(1))
    if (anyDuplicated(observed_names) > 0L ||
        !identical(sort(observed_names), expected_names)) {
        protDiaEnrichArtifactAbort(
            "enrichment product set is incomplete or unexpected",
            "multischolar_prot_dia_enrichment_product_mismatch"
        )
    }
    products
}

protDiaEnrichValidateResultSemantics <- function(table, parameters, records) {
    required <- c("directionality", "analysis_method")
    if (!is.data.frame(table) || !all(required %in% names(table))) {
        protDiaEnrichArtifactAbort(
            "enrichment display result schema is incomplete",
            "multischolar_invalid_prot_dia_enrichment_table"
        )
    }
    expected_method <- parameters$backend
    expected_labels <- if (identical(expected_method, "gprofiler2")) {
        c(up = "positive", down = "negative")
    } else {
        c(up = "up", down = "down")
    }
    expected_rows <- stats::setNames(integer(2L), names(expected_labels))
    for (record in records) {
        if (identical(record$status, "succeeded")) {
            expected_rows[[record$direction]] <- record$response$rows
        }
    }
    observed_rows <- vapply(expected_labels, function(label) {
        sum(table$directionality == label, na.rm = TRUE)
    }, integer(1))
    valid <- all(table$analysis_method == expected_method) &&
        all(table$directionality %in% unname(expected_labels)) &&
        identical(observed_rows, expected_rows) &&
        identical(nrow(table), sum(expected_rows))
    if (!isTRUE(valid)) {
        protDiaEnrichArtifactAbort(
            "enrichment display results differ from service provenance",
            "multischolar_prot_dia_enrichment_result_mismatch"
        )
    }
    invisible(table)
}
