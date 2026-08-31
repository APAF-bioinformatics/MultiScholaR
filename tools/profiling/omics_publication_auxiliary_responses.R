auxPublicationResponseDefinitions <- function() {
    result_fields <- list(
        "termDescription", "enrichmentScore", "falseDiscoveryRate",
        "genesMapped", "mappedIDs", "comparison"
    )
    define <- function(service_id, request_schema, response_schema) {
        request <- list(
            request_id = paste0(service_id, ".synthetic.request.v1"),
            schema = request_schema,
            service_version = "contract_fixture_v1",
            identifier_namespace = if (identical(service_id, "stringdb")) {
                "synthetic_protein_id"
            } else {
                "synthetic_chebi_id"
            }
        )
        response <- list(
            response_id = paste0(service_id, ".synthetic.response.v1"),
            schema = response_schema,
            service_version = "contract_fixture_v1",
            result_fields = result_fields,
            row_count = 5L,
            content_origin = "repository_generated_schema_fixture"
        )
        list(
            service_id = service_id,
            request = request,
            response = response,
            request_response_sha256 = publicationObjectDigest(list(
                request = request,
                response = response
            )),
            redistribution_status = "self_generated_no_external_content",
            service_behaviour_authority = FALSE,
            biological_authority = FALSE,
            local_transform_performance_authority = TRUE,
            live_probe_allowed = FALSE,
            promotion_authority = FALSE
        )
    }
    list(
        define(
            "stringdb",
            "multischolar.synthetic_string_rank_request/1.0.0",
            "multischolar.synthetic_string_enrichment_response/1.0.0"
        ),
        define(
            "kegg",
            "multischolar.synthetic_kegg_rank_request/1.0.0",
            "multischolar.synthetic_kegg_enrichment_response/1.0.0"
        ),
        define(
            "reactome",
            "multischolar.synthetic_reactome_rank_request/1.0.0",
            "multischolar.synthetic_reactome_enrichment_response/1.0.0"
        ),
        define(
            "metabolomics_pathway",
            "multischolar.synthetic_metabolite_rank_request/1.0.0",
            "multischolar.synthetic_metabolite_enrichment_response/1.0.0"
        )
    )
}

auxPublicationValidateResponseDefinition <- function(definition) {
    publicationRequireNames(definition, c(
        "service_id", "request", "response", "request_response_sha256",
        "redistribution_status", "service_behaviour_authority",
        "biological_authority", "local_transform_performance_authority",
        "live_probe_allowed", "promotion_authority"
    ), "Auxiliary response definition")
    publicationRequireNames(definition$request, c(
        "request_id", "schema", "service_version", "identifier_namespace"
    ), "Auxiliary response request")
    publicationRequireNames(definition$response, c(
        "response_id", "schema", "service_version", "result_fields",
        "row_count", "content_origin"
    ), "Auxiliary response payload")
    expected <- Filter(function(value) {
        identical(value$service_id, definition$service_id)
    }, auxPublicationResponseDefinitions())
    digest <- publicationObjectDigest(list(
        request = definition$request,
        response = definition$response
    ))
    valid <- length(expected) == 1L && identical(definition, expected[[1L]]) &&
        identical(definition$request_response_sha256, digest) &&
        identical(
            definition$redistribution_status,
            "self_generated_no_external_content"
        ) &&
        !isTRUE(definition$service_behaviour_authority) &&
        !isTRUE(definition$biological_authority) &&
        isTRUE(definition$local_transform_performance_authority) &&
        !isTRUE(definition$live_probe_allowed) &&
        !isTRUE(definition$promotion_authority)
    if (!valid) auxPublicationAbort("auxiliary response definition differs")
    invisible(definition)
}

auxPublicationBuildResponseAuthority <- function() {
    list(
        schema = "multischolar.omics_publication_auxiliary_responses",
        schema_version = .AUX_PUBLICATION_VERSION,
        responses_id = paste0(
            "multischolar.omics_publication_auxiliary_responses.2026-08-28.v1"
        ),
        owner_ticket_id = .AUX_PUBLICATION_OWNER,
        status = "frozen_pre_candidate",
        responses = auxPublicationResponseDefinitions(),
        primary_network_allowed = FALSE,
        live_checks_separate = TRUE,
        live_checks_promotional = FALSE,
        generic_cross_service_fixture_allowed = FALSE,
        unknown_policy = "reject",
        publication_authority = FALSE
    )
}

auxPublicationValidateResponses <- function(record) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "responses_id", "owner_ticket_id", "status",
        "responses", "primary_network_allowed", "live_checks_separate",
        "live_checks_promotional", "generic_cross_service_fixture_allowed",
        "unknown_policy", "publication_authority"
    ), "Auxiliary response authority")
    lapply(record$responses, auxPublicationValidateResponseDefinition)
    ids <- vapply(record$responses, `[[`, character(1), "service_id")
    expected_ids <- vapply(
        auxPublicationResponseDefinitions(),
        `[[`,
        character(1),
        "service_id"
    )
    valid <- identical(
        record$schema,
        "multischolar.omics_publication_auxiliary_responses"
    ) && identical(record$schema_version, .AUX_PUBLICATION_VERSION) &&
        identical(record$owner_ticket_id, .AUX_PUBLICATION_OWNER) &&
        identical(record$status, "frozen_pre_candidate") &&
        identical(ids, expected_ids) && !anyDuplicated(ids) &&
        !isTRUE(record$primary_network_allowed) &&
        isTRUE(record$live_checks_separate) &&
        !isTRUE(record$live_checks_promotional) &&
        !isTRUE(record$generic_cross_service_fixture_allowed) &&
        identical(record$unknown_policy, "reject") &&
        !isTRUE(record$publication_authority)
    if (!valid) auxPublicationAbort("auxiliary response authority differs")
    invisible(record)
}

auxPublicationResponseDefinition <- function(responses, service_id) {
    matches <- Filter(function(value) {
        identical(value$service_id, service_id)
    }, responses$responses)
    if (length(matches) != 1L) {
        auxPublicationAbort(
            "auxiliary response service is missing or duplicated",
            "response_excluded"
        )
    }
    matches[[1L]]
}

auxPublicationExpandResponse <- function(
    definition,
    row_count,
    identifier_multiplicity,
    pathway_count
) {
    auxPublicationValidateResponseDefinition(definition)
    index <- seq_len(row_count)
    mapped <- vapply(index, function(value) {
        ids <- seq.int(
            (value - 1L) * identifier_multiplicity + 1L,
            value * identifier_multiplicity
        )
        if (identical(definition$service_id, "stringdb")) {
            return(paste0(
                "SYNPROT",
                sprintf("%09d", ids),
                collapse = ","
            ))
        }
        paste0("CHEBI:", 1000000L + ids, collapse = ",")
    }, character(1))
    data.frame(
        termDescription = paste(
            "Synthetic",
            definition$service_id,
            "pathway",
            sprintf("%09d", ((index - 1L) %% pathway_count) + 1L)
        ),
        enrichmentScore = round(log1p(index), digits = 8L),
        falseDiscoveryRate = pmin(1, index / (row_count * 20)),
        genesMapped = rep.int(as.integer(identifier_multiplicity), row_count),
        mappedIDs = mapped,
        comparison = paste0(definition$service_id, "_local_fixture"),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
}

auxPublicationValidateExpandedResponse <- function(
    response,
    definition,
    row_count,
    identifier_multiplicity,
    pathway_count
) {
    expected <- auxPublicationExpandResponse(
        definition,
        row_count,
        identifier_multiplicity,
        pathway_count
    )
    valid <- is.data.frame(response) && identical(names(response), names(expected)) &&
        identical(nrow(response), as.integer(row_count)) &&
        identical(
            publicationObjectDigest(response),
            publicationObjectDigest(expected)
        ) && !anyNA(response$mappedIDs) && all(nzchar(response$mappedIDs))
    if (!valid) {
        auxPublicationAbort(
            "auxiliary expanded response differs",
            "response_invalid"
        )
    }
    invisible(response)
}
