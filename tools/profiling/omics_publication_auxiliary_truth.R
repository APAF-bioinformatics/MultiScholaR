auxPublicationContractBasis <- function(contract) {
    basis <- contract
    basis$expected_digests <- list(
        payload_sha256 = strrep("0", 64L),
        truth_sha256 = strrep("0", 64L)
    )
    publicationObjectDigest(basis)
}

auxPublicationCountValues <- function(values, levels) {
    stats::setNames(lapply(levels, function(level) {
        as.integer(sum(values == level, na.rm = TRUE))
    }), levels)
}

auxPublicationPhosphositeTruth <- function(payload, dimensions) {
    rows <- payload$row_plan
    valid_accession <- rows$category_id %in% c("single", "ambiguous")
    accession_multiplicity <- ifelse(rows$category_id == "ambiguous", 2L, 1L)
    abundance_eligible <- !rows$zero_abundance
    site_probability_rows <- sum(
        accession_multiplicity[abundance_eligible & valid_accession]
    )
    mapping_rows <- rep(seq_len(nrow(rows)), accession_multiplicity)
    mapping_proteins <- unlist(lapply(seq_len(nrow(rows)), function(index) {
        if (identical(rows$category_id[[index]], "ambiguous")) {
            c(rows$protein_index[[index]], rows$second_protein_index[[index]])
        } else {
            rows$protein_index[[index]]
        }
    }), use.names = FALSE)
    mapping_valid <- abundance_eligible[mapping_rows] &
        valid_accession[mapping_rows]
    high <- mapping_valid & rows$high_confidence[mapping_rows]
    high_any <- tapply(high, mapping_proteins, any)
    high_multisite <- tapply(
        high & rows$multisite[mapping_rows],
        mapping_proteins,
        any
    )
    protein_key <- as.character(mapping_proteins)
    rescued <- ifelse(
        rows$multisite[mapping_rows],
        high_multisite[protein_key],
        high_any[protein_key]
    )
    rescued[is.na(rescued)] <- FALSE
    filtered_rows <- sum(mapping_valid & rescued)
    filtered_evidence_rows <- length(unique(mapping_rows[
        mapping_valid & rescued
    ]))
    list(
        evidence_row_count = as.integer(nrow(payload$evidence)),
        protein_count = as.integer(nrow(payload$proteins)),
        accession_category_counts = auxPublicationCountValues(
            rows$category_id,
            c("single", "ambiguous", "missing", "contaminant", "reverse")
        ),
        duplicate_evidence_count = as.integer(sum(rows$duplicate)),
        multisite_evidence_count = as.integer(sum(rows$multisite)),
        high_confidence_evidence_count = as.integer(sum(rows$high_confidence)),
        zero_abundance_evidence_count = as.integer(sum(rows$zero_abundance)),
        repeated_peptide_protein_count = as.integer(sum(
            payload$proteins$repeated_peptide
        )),
        expected_accession_rows = as.integer(sum(accession_multiplicity)),
        expected_after_abundance_filter = as.integer(sum(abundance_eligible)),
        expected_site_probability_rows = as.integer(site_probability_rows),
        expected_score_filtered_rows = as.integer(filtered_rows),
        expected_long_rows = as.numeric(
            filtered_rows * dimensions$measured_sample_count
        ),
        expected_paralog_rows = as.numeric(
            filtered_evidence_rows * dimensions$measured_sample_count
        ),
        measured_sample_count = as.integer(dimensions$measured_sample_count),
        input_column_order = as.list(names(payload$evidence)),
        known_deprecations = list("purrr.cross2"),
        exact_order_required = TRUE
    )
}

auxPublicationLayerJoinCardinality <- function(weights, assay_data, view) {
    selected <- weights[weights$view == view, , drop = FALSE]
    selected$feature <- sub(
        paste0("_", view, "$"),
        "",
        selected$feature
    )
    assay <- assay_data[[view]]
    weight_counts <- table(selected$feature, useNA = "no")
    assay_counts <- table(assay$metabolite, useNA = "no")
    shared <- intersect(names(weight_counts), names(assay_counts))
    as.numeric(sum(weight_counts[shared] * assay_counts[shared]))
}

auxPublicationTopWeightIds <- function(matrix, limit = 20L) {
    values <- matrix[, 1L]
    ordered <- order(-abs(values), na.last = TRUE, method = "radix")
    normalized_row_names <- rownames(as.data.frame(matrix))
    normalized_row_names[utils::head(ordered, limit)]
}

auxPublicationMultiomicsTruth <- function(payload, dimensions) {
    feature_plan <- payload$feature_plan
    layer_names <- c("proteome", "metabolome_lc", "metabolome_gc")
    layer_counts <- stats::setNames(lapply(layer_names, function(layer) {
        as.integer(sum(feature_plan$layer == layer))
    }), layer_names)
    sample_union <- sort(unique(unlist(
        payload$sample_registry,
        use.names = FALSE
    )), method = "radix")
    sample_intersection <- Reduce(intersect, payload$sample_registry)
    top_ids <- if (length(payload$model_double$weights)) {
        stats::setNames(lapply(
            payload$model_double$weights,
            auxPublicationTopWeightIds
        ), layer_names)
    } else {
        list()
    }
    join_cardinality <- if (length(payload$assay_data)) {
        stats::setNames(lapply(
            c("metabolome_lc", "metabolome_gc"),
            function(view) {
                auxPublicationLayerJoinCardinality(
                    payload$weights,
                    payload$assay_data,
                    view
                )
            }
        ), c("metabolome_lc", "metabolome_gc"))
    } else {
        list()
    }
    response_counts <- if (length(payload$responses)) {
        lapply(payload$responses, nrow)
    } else {
        list()
    }
    response_digests <- if (length(payload$responses)) {
        lapply(payload$responses, publicationObjectDigest)
    } else {
        list()
    }
    list(
        feature_count = as.integer(nrow(feature_plan)),
        layer_feature_counts = layer_counts,
        missing_weight_count = as.integer(sum(feature_plan$missing)),
        duplicate_feature_count = as.integer(sum(feature_plan$duplicated)),
        tied_weight_count = as.integer(sum(feature_plan$tied)),
        contextual_sample_count = as.integer(dimensions$sample_count),
        measured_sample_count = as.integer(dimensions$measured_sample_count),
        sample_union_count = as.integer(length(sample_union)),
        sample_intersection_count = as.integer(length(sample_intersection)),
        sample_union_sha256 = publicationObjectDigest(as.list(sample_union)),
        sample_intersection_sha256 = publicationObjectDigest(
            as.list(sort(sample_intersection, method = "radix"))
        ),
        top_weight_ids = top_ids,
        top_weight_ids_sha256 = lapply(top_ids, publicationObjectDigest),
        enrichment_join_cardinality = join_cardinality,
        response_row_counts = response_counts,
        response_sha256 = response_digests,
        sample_scaling_performance_claim = FALSE,
        mofa_fitting_or_inference_claim = FALSE,
        service_behaviour_claim = FALSE,
        exact_order_required = TRUE
    )
}

auxPublicationBuildTruth <- function(contract, payload, payload_sha256) {
    auxPublicationValidateContract(contract)
    route_id <- contract$route$route_id
    facts <- if (identical(route_id, "phosphosite_stages")) {
        auxPublicationPhosphositeTruth(payload, contract$dimensions)
    } else {
        auxPublicationMultiomicsTruth(payload, contract$dimensions)
    }
    list(
        schema = "multischolar.omics_publication_auxiliary_truth",
        schema_version = .AUX_PUBLICATION_VERSION,
        truth_id = paste0(contract$workload_id, ".truth"),
        owner_ticket_id = .AUX_PUBLICATION_OWNER,
        status = "frozen_pre_candidate",
        workload_id = contract$workload_id,
        route_id = route_id,
        workload_class = contract$workload_class,
        contract_basis_sha256 = auxPublicationContractBasis(contract),
        payload_sha256 = payload_sha256,
        dimensions = contract$dimensions,
        facts = facts,
        oracle = list(
            implementation = "independent_generator_cardinality_and_order_oracle",
            package_reader_used = FALSE,
            package_export_used = FALSE,
            historical_result_used = FALSE,
            candidate_result_used = FALSE
        ),
        scientific_authority = FALSE,
        promotion_authority = FALSE,
        publication_authority = FALSE
    )
}

auxPublicationValidateTruth <- function(truth, contract, payload_path) {
    publicationRequireNames(truth, c(
        "schema", "schema_version", "truth_id", "owner_ticket_id", "status",
        "workload_id", "route_id", "workload_class", "contract_basis_sha256",
        "payload_sha256", "dimensions", "facts", "oracle",
        "scientific_authority", "promotion_authority", "publication_authority"
    ), "Auxiliary truth")
    publicationRequireNames(truth$oracle, c(
        "implementation", "package_reader_used", "package_export_used",
        "historical_result_used", "candidate_result_used"
    ), "Auxiliary truth oracle")
    auxPublicationValidateContract(contract)
    valid <- identical(
        truth$schema,
        "multischolar.omics_publication_auxiliary_truth"
    ) && identical(truth$schema_version, .AUX_PUBLICATION_VERSION) &&
        identical(truth$owner_ticket_id, .AUX_PUBLICATION_OWNER) &&
        identical(truth$status, "frozen_pre_candidate") &&
        identical(truth$workload_id, contract$workload_id) &&
        identical(truth$route_id, contract$route$route_id) &&
        identical(truth$workload_class, contract$workload_class) &&
        identical(truth$contract_basis_sha256, auxPublicationContractBasis(contract)) &&
        identical(truth$payload_sha256, publicationFileDigest(payload_path)) &&
        identical(truth$payload_sha256, contract$expected_digests$payload_sha256) &&
        identical(
            publicationObjectDigest(truth$dimensions),
            publicationObjectDigest(contract$dimensions)
        ) &&
        !isTRUE(truth$oracle$package_reader_used) &&
        !isTRUE(truth$oracle$package_export_used) &&
        !isTRUE(truth$oracle$historical_result_used) &&
        !isTRUE(truth$oracle$candidate_result_used) &&
        !isTRUE(truth$scientific_authority) &&
        !isTRUE(truth$promotion_authority) &&
        !isTRUE(truth$publication_authority)
    if (!valid) auxPublicationAbort("auxiliary truth differs")
    invisible(truth)
}

auxPublicationRecalculateTruth <- function(contract, payload_path) {
    payload <- readRDS(payload_path)
    auxPublicationBuildTruth(
        contract,
        payload,
        publicationFileDigest(payload_path)
    )
}

auxPublicationExpectedTopIds <- function(truth, view) {
    ids <- unlist(truth$facts$top_weight_ids[[view]], use.names = FALSE)
    if (identical(view, "proteome")) {
        return(sub("_proteome$", "", ids))
    }
    if (identical(view, "transcriptome")) {
        return(sub("_transcriptome$", "", ids))
    }
    ids
}

auxPublicationValidatePhosphositeResult <- function(result, truth) {
    expected <- truth$facts
    stages <- result$stages
    valid <- identical(stages$cleaned, expected$evidence_row_count) &&
        identical(stages$accession, expected$expected_accession_rows) &&
        identical(stages$abundance, expected$expected_after_abundance_filter) &&
        identical(
            stages$probabilities,
            expected$expected_site_probability_rows
        ) &&
        identical(stages$filtered, expected$expected_score_filtered_rows) &&
        identical(
            as.numeric(stages$long),
            as.numeric(expected$expected_long_rows)
        ) && identical(
            as.numeric(stages$paralog),
            as.numeric(expected$expected_paralog_rows)
        )
    if (!valid) auxPublicationAbort("phosphosite stage result differs")
    invisible(result)
}

auxPublicationValidateMofaResult <- function(result, truth) {
    plots <- result$retained
    views <- names(truth$facts$layer_feature_counts)
    plot_ids <- lapply(views, function(view) {
        as.character(plots[[view]]$data$gene_symbol)
    })
    names(plot_ids) <- views
    valid <- identical(result$stages$view_count, as.integer(length(views))) &&
        identical(
            result$stages$ranked_feature_count,
            truth$facts$feature_count
        ) &&
        identical(
            result$stages$plotted_feature_count,
            as.integer(length(views) * 20L)
        ) &&
        all(vapply(views, function(view) {
            identical(plot_ids[[view]], auxPublicationExpectedTopIds(truth, view))
        }, logical(1)))
    if (!valid) auxPublicationAbort("MOFA weight result differs")
    invisible(result)
}

auxPublicationValidateEnrichmentResult <- function(result, truth) {
    joins <- unlist(
        truth$facts$enrichment_join_cardinality,
        use.names = FALSE
    )
    expected_lengths <- as.list(rep(joins, each = 2L))
    response_rows <- unlist(truth$facts$response_row_counts, use.names = FALSE)
    expected_rows <- as.integer(2L * sum(response_rows))
    valid <- identical(
        result$stages$response_service_order,
        list("kegg", "reactome", "kegg", "reactome")
    ) && identical(
        as.numeric(unlist(result$stages$ranked_list_lengths, use.names = FALSE)),
        as.numeric(unlist(expected_lengths, use.names = FALSE))
    ) &&
        identical(result$stages$result_row_count, expected_rows)
    if (!valid) auxPublicationAbort("metabolite enrichment result differs")
    invisible(result)
}

auxPublicationValidatePathwayResult <- function(result, truth) {
    response_rows <- unlist(truth$facts$response_row_counts, use.names = FALSE)
    expected_rows <- as.integer(2L * sum(response_rows))
    valid <- identical(
        result$stages$assay_call_order,
        list("metabolome_lc", "metabolome_gc")
    ) && identical(result$stages$result_row_count, expected_rows) &&
        identical(result$stages$mapped_name_count, expected_rows)
    if (!valid) auxPublicationAbort("metabolite pathway result differs")
    invisible(result)
}

auxPublicationValidateStringResult <- function(result, truth) {
    response_rows <- unlist(truth$facts$response_row_counts, use.names = FALSE)
    valid <- identical(
        result$stages$ranked_identifier_count,
        truth$facts$feature_count
    ) && identical(
        result$stages$result_row_count,
        as.integer(sum(response_rows))
    )
    if (!valid) auxPublicationAbort("STRING rank result differs")
    invisible(result)
}

auxPublicationValidateRunResult <- function(result, truth, contract) {
    if (!identical(truth$route_id, contract$route$route_id)) {
        auxPublicationAbort("auxiliary result route differs")
    }
    switch(
        contract$route$route_id,
        phosphosite_stages = auxPublicationValidatePhosphositeResult(
            result,
            truth
        ),
        mofa_weights = auxPublicationValidateMofaResult(result, truth),
        metabolite_enrichment = auxPublicationValidateEnrichmentResult(
            result,
            truth
        ),
        metabolite_pathway = auxPublicationValidatePathwayResult(result, truth),
        stringdb_rank = auxPublicationValidateStringResult(result, truth),
        auxPublicationAbort("auxiliary result validator route is unsupported")
    )
    invisible(result)
}
