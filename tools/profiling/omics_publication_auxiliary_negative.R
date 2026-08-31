auxPublicationNegativeDefinitions <- function() {
    define <- function(case_id, target, mutation, expected_class) {
        list(
            case_id = case_id,
            target = target,
            mutation = mutation,
            expected_condition_class = expected_class,
            expected_rejection_stage = "admission",
            network_invocation_allowed = FALSE,
            runner_invocation_allowed = FALSE,
            artifact_invocation_allowed = FALSE,
            backend_invocation_allowed = FALSE,
            promotion_authority = FALSE
        )
    }
    contract_error <- "multischolar_publication_auxiliary_contract_error"
    response_error <- "multischolar_publication_auxiliary_response_invalid"
    list(
        define("contract.unknown_route", "contract", "unknown_route", contract_error),
        define(
            "contract.network_enabled",
            "contract",
            "network_enabled",
            contract_error
        ),
        define(
            "contract.workflow_authority",
            "contract",
            "workflow_authority",
            contract_error
        ),
        define(
            "contract.artifact_authority",
            "contract",
            "artifact_authority",
            contract_error
        ),
        define(
            "contract.mofa_fitting_authority",
            "contract",
            "mofa_fitting_authority",
            contract_error
        ),
        define(
            "contract.stale_parameter_binding",
            "contract",
            "stale_parameter_binding",
            contract_error
        ),
        define(
            "definition.wrong_service",
            "response_definition",
            "wrong_service",
            contract_error
        ),
        define(
            "definition.wrong_version",
            "response_definition",
            "wrong_version",
            contract_error
        ),
        define(
            "definition.stale_digest",
            "response_definition",
            "stale_digest",
            contract_error
        ),
        define(
            "authority.generic_cross_service",
            "response_authority",
            "generic_cross_service",
            contract_error
        ),
        define(
            "response.reordered",
            "expanded_response",
            "reordered",
            response_error
        ),
        define(
            "response.missing_ids",
            "expanded_response",
            "missing_ids",
            response_error
        ),
        define(
            "response.truncated",
            "expanded_response",
            "truncated",
            response_error
        ),
        define(
            "response.schema_drift",
            "expanded_response",
            "schema_drift",
            response_error
        ),
        define(
            "response.duplicate_row",
            "expanded_response",
            "duplicate_row",
            response_error
        ),
        define(
            "split.pilot_source_overlap",
            "split",
            "pilot_source_overlap",
            contract_error
        ),
        define(
            "split.candidate_access",
            "split",
            "candidate_access",
            contract_error
        ),
        define(
            "phosphosite.malformed_probability",
            "phosphosite_payload",
            "malformed_probability",
            contract_error
        ),
        define(
            "phosphosite.missing_sequence",
            "phosphosite_payload",
            "missing_sequence",
            contract_error
        ),
        define(
            "multiomics.unknown_layer",
            "multiomics_payload",
            "unknown_layer",
            contract_error
        ),
        define(
            "multiomics.wrong_response_service",
            "multiomics_payload",
            "wrong_response_service",
            contract_error
        ),
        define(
            "service.unknown",
            "response_lookup",
            "unknown_service",
            "multischolar_publication_auxiliary_response_excluded"
        )
    )
}

auxPublicationBuildNegativeAuthority <- function() {
    list(
        schema = "multischolar.omics_publication_auxiliary_negatives",
        schema_version = .AUX_PUBLICATION_VERSION,
        negatives_id = paste0(
            "multischolar.omics_publication_auxiliary_negatives.2026-08-28.v1"
        ),
        owner_ticket_id = .AUX_PUBLICATION_OWNER,
        status = "frozen_pre_candidate",
        surface_binding = auxPublicationBinding(paste0(
            "tests/testdata/omics-performance/auxiliary/surfaces-v1.json"
        )),
        exclusion_binding = auxPublicationBinding(paste0(
            "tests/testdata/omics-performance/auxiliary/exclusions-v1.json"
        )),
        response_binding = auxPublicationBinding(paste0(
            "tests/testdata/omics-performance/auxiliary/responses-v1.json"
        )),
        split_binding = auxPublicationBinding(paste0(
            "tests/testdata/omics-performance/auxiliary/splits-v1.json"
        )),
        cases = auxPublicationNegativeDefinitions(),
        unknown_policy = "reject",
        publication_authority = FALSE
    )
}

auxPublicationValidateNegativeCase <- function(case) {
    publicationRequireNames(case, c(
        "case_id", "target", "mutation", "expected_condition_class",
        "expected_rejection_stage", "network_invocation_allowed",
        "runner_invocation_allowed", "artifact_invocation_allowed",
        "backend_invocation_allowed", "promotion_authority"
    ), "Auxiliary negative case")
    expected <- Filter(function(value) {
        identical(value$case_id, case$case_id)
    }, auxPublicationNegativeDefinitions())
    valid <- length(expected) == 1L && identical(case, expected[[1L]]) &&
        identical(case$expected_rejection_stage, "admission") &&
        !isTRUE(case$network_invocation_allowed) &&
        !isTRUE(case$runner_invocation_allowed) &&
        !isTRUE(case$artifact_invocation_allowed) &&
        !isTRUE(case$backend_invocation_allowed) &&
        !isTRUE(case$promotion_authority)
    if (!valid) auxPublicationAbort("auxiliary negative case differs")
    invisible(case)
}

auxPublicationValidateNegatives <- function(record) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "negatives_id", "owner_ticket_id", "status",
        "surface_binding", "exclusion_binding", "response_binding",
        "split_binding", "cases", "unknown_policy", "publication_authority"
    ), "Auxiliary negative authority")
    lapply(c(
        "surface_binding", "exclusion_binding", "response_binding",
        "split_binding"
    ), function(name) {
        auxPublicationValidateBinding(record[[name]], name)
    })
    lapply(record$cases, auxPublicationValidateNegativeCase)
    ids <- vapply(record$cases, `[[`, character(1), "case_id")
    expected_ids <- vapply(
        auxPublicationNegativeDefinitions(),
        `[[`,
        character(1),
        "case_id"
    )
    valid <- identical(
        record$schema,
        "multischolar.omics_publication_auxiliary_negatives"
    ) && identical(record$schema_version, .AUX_PUBLICATION_VERSION) &&
        identical(record$owner_ticket_id, .AUX_PUBLICATION_OWNER) &&
        identical(record$status, "frozen_pre_candidate") &&
        identical(ids, expected_ids) && !anyDuplicated(ids) &&
        identical(record$unknown_policy, "reject") &&
        !isTRUE(record$publication_authority)
    if (!valid) auxPublicationAbort("auxiliary negative authority differs")
    invisible(record)
}

auxPublicationNegativeContract <- function(route_id) {
    splits <- publicationReadJson(paste0(
        "tests/testdata/omics-performance/auxiliary/splits-v1.json"
    ))
    assignment <- auxPublicationFindAssignment(
        splits,
        route_id,
        "fixture_correctness"
    )
    auxPublicationBuildContract(
        route_id,
        "fixture_correctness",
        assignment,
        NULL
    )
}

auxPublicationEvaluateContractNegative <- function(mutation) {
    contract <- auxPublicationNegativeContract("mofa_weights")
    if (identical(mutation, "unknown_route")) {
        contract$route$route_id <- "unknown.route"
    } else if (identical(mutation, "network_enabled")) {
        contract$execution$network_allowed <- TRUE
    } else if (identical(mutation, "workflow_authority")) {
        contract$claim_scope$workflow_authority <- TRUE
    } else if (identical(mutation, "artifact_authority")) {
        contract$claim_scope$artifact_authority <- TRUE
    } else if (identical(mutation, "mofa_fitting_authority")) {
        contract$truth_contract$mofa_fitting_authority <- TRUE
    } else if (identical(mutation, "stale_parameter_binding")) {
        contract$parameter_authority$sha256 <- strrep("0", 64L)
    }
    auxPublicationValidateContract(contract)
}

auxPublicationEvaluateDefinitionNegative <- function(mutation) {
    responses <- publicationReadJson(paste0(
        "tests/testdata/omics-performance/auxiliary/responses-v1.json"
    ))
    definition <- auxPublicationResponseDefinition(responses, "kegg")
    if (identical(mutation, "wrong_service")) {
        definition$service_id <- "wrong_service"
    } else if (identical(mutation, "wrong_version")) {
        definition$request$service_version <- "wrong_version"
    } else {
        definition$request_response_sha256 <- strrep("0", 64L)
    }
    auxPublicationValidateResponseDefinition(definition)
}

auxPublicationEvaluateExpandedResponseNegative <- function(mutation) {
    responses <- publicationReadJson(paste0(
        "tests/testdata/omics-performance/auxiliary/responses-v1.json"
    ))
    definition <- auxPublicationResponseDefinition(responses, "kegg")
    response <- auxPublicationExpandResponse(definition, 12L, 4L, 250L)
    if (identical(mutation, "reordered")) {
        response <- response[rev(seq_len(nrow(response))), , drop = FALSE]
    } else if (identical(mutation, "missing_ids")) {
        response$mappedIDs[[1L]] <- NA_character_
    } else if (identical(mutation, "truncated")) {
        response <- response[-nrow(response), , drop = FALSE]
    } else if (identical(mutation, "schema_drift")) {
        names(response)[[1L]] <- "term_name"
    } else {
        response[2L, ] <- response[1L, ]
    }
    auxPublicationValidateExpandedResponse(
        response,
        definition,
        12L,
        4L,
        250L
    )
}

auxPublicationEvaluateSplitNegative <- function(mutation) {
    splits <- publicationReadJson(paste0(
        "tests/testdata/omics-performance/auxiliary/splits-v1.json"
    ))
    if (identical(mutation, "pilot_source_overlap")) {
        splits$pilot_calibration_assignments[[1L]]$source_id <-
            splits$assignments[[1L]]$source_id
    } else {
        splits$candidate_access_allowed <- TRUE
    }
    auxPublicationValidateSplits(splits)
}

auxPublicationEvaluatePayloadNegative <- function(target, mutation) {
    route_id <- if (identical(target, "phosphosite_payload")) {
        "phosphosite_stages"
    } else {
        "metabolite_enrichment"
    }
    contract <- auxPublicationNegativeContract(route_id)
    payload <- auxPublicationBuildPayload(contract)
    if (identical(mutation, "malformed_probability")) {
        payload$evidence$phospho_sty_probabilities[[1L]] <- "MALFORMED"
    } else if (identical(mutation, "missing_sequence")) {
        payload$proteins$seq[[1L]] <- NA_character_
    } else if (identical(mutation, "unknown_layer")) {
        payload$feature_plan$layer[[1L]] <- "unknown_layer"
    } else {
        names(payload$responses)[[1L]] <- "wrong_service"
    }
    auxPublicationValidatePayload(payload, contract)
}

auxPublicationEvaluateNegative <- function(case) {
    auxPublicationValidateNegativeCase(case)
    error <- tryCatch({
        switch(
            case$target,
            contract = auxPublicationEvaluateContractNegative(case$mutation),
            response_definition = auxPublicationEvaluateDefinitionNegative(
                case$mutation
            ),
            response_authority = {
                record <- publicationReadJson(paste0(
                    "tests/testdata/omics-performance/auxiliary/",
                    "responses-v1.json"
                ))
                record$generic_cross_service_fixture_allowed <- TRUE
                auxPublicationValidateResponses(record)
            },
            expanded_response = auxPublicationEvaluateExpandedResponseNegative(
                case$mutation
            ),
            split = auxPublicationEvaluateSplitNegative(case$mutation),
            phosphosite_payload = auxPublicationEvaluatePayloadNegative(
                case$target,
                case$mutation
            ),
            multiomics_payload = auxPublicationEvaluatePayloadNegative(
                case$target,
                case$mutation
            ),
            response_lookup = {
                responses <- publicationReadJson(paste0(
                    "tests/testdata/omics-performance/auxiliary/",
                    "responses-v1.json"
                ))
                auxPublicationResponseDefinition(responses, "unknown_service")
            }
        )
        NULL
    }, error = function(error) error)
    if (!inherits(error, case$expected_condition_class)) {
        auxPublicationAbort("auxiliary negative outcome differs")
    }
    list(
        case_id = case$case_id,
        condition_class = class(error)[[1L]],
        network_invocations = 0L,
        runner_invocations = 0L,
        artifact_invocations = 0L,
        backend_invocations = 0L,
        promotion_authority = FALSE
    )
}
