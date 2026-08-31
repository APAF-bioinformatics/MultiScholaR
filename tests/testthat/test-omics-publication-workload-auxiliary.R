.publicationAuxRepoPath <- function(...) {
    file.path(
        normalizePath(testthat::test_path("..", ".."), mustWork = TRUE),
        ...
    )
}

for (.publication_aux_source in c(
    "omics_publication_protocol.R",
    "omics_publication_auxiliary_contracts.R",
    "omics_publication_auxiliary_model.R",
    "omics_publication_auxiliary_responses.R",
    "omics_publication_auxiliary_sources.R",
    "omics_publication_auxiliary_governance.R",
    "omics_publication_auxiliary_truth.R",
    "omics_publication_workload_auxiliary.R",
    "omics_publication_auxiliary_negative.R",
    "omics_publication_auxiliary_manifest.R"
)) {
    sys.source(
        .publicationAuxRepoPath(
            "tools",
            "profiling",
            .publication_aux_source
        ),
        envir = environment()
    )
}
rm(.publication_aux_source)

.publicationAuxRecord <- function(name) {
    publicationReadJson(.publicationAuxRepoPath(
        "tests",
        "testdata",
        "omics-performance",
        "auxiliary",
        name
    ))
}

test_that("auxiliary surfaces exactly bind current exported APIs", {
    authority <- .publicationAuxRecord("surfaces-v1.json")

    expect_silent(auxPublicationValidateSurfaceAuthority(authority))
    expect_identical(
        publicationObjectDigest(authority),
        publicationObjectDigest(auxPublicationBuildSurfaceAuthority())
    )
    expect_identical(
        vapply(authority$surfaces, `[[`, character(1), "surface_id"),
        c("phosphosite.api", "multiomics.api")
    )

    entry_points <- unlist(lapply(authority$surfaces, function(surface) {
        c(surface$primary_entry_point, surface$compatibility_entry_points)
    }), use.names = FALSE)
    expect_identical(anyDuplicated(entry_points), 0L)
    expect_true(all(entry_points %in% getNamespaceExports("MultiScholaR")))
    expect_true(all(vapply(entry_points, function(name) {
        is.function(getExportedValue("MultiScholaR", name))
    }, logical(1))))
})

test_that("auxiliary surfaces cannot acquire workflow or backend authority", {
    authority <- .publicationAuxRecord("surfaces-v1.json")

    workflow <- publicationGovernanceCopy(authority)
    workflow$workflow_capability_creation_allowed <- TRUE
    expect_error(
        auxPublicationValidateSurfaceAuthority(workflow),
        class = "multischolar_publication_auxiliary_contract_error"
    )

    backend <- publicationGovernanceCopy(authority)
    backend$automatic_backend_decision_allowed <- TRUE
    expect_error(
        auxPublicationValidateSurfaceAuthority(backend),
        class = "multischolar_publication_auxiliary_contract_error"
    )

    widened <- publicationGovernanceCopy(authority)
    widened$surfaces[[1L]]$claim_scope <- "publication_and_promotion"
    expect_error(
        auxPublicationValidateSurfaceAuthority(widened),
        class = "multischolar_publication_auxiliary_contract_error"
    )

    stale <- publicationGovernanceCopy(authority)
    stale$capability_binding$sha256 <- strrep("0", 64L)
    expect_error(
        auxPublicationValidateSurfaceAuthority(stale),
        class = "multischolar_publication_auxiliary_contract_error"
    )
})

test_that("auxiliary surface resolution rejects unknown identities first", {
    calls <- integer(4L)
    names(calls) <- c("prepare", "runner", "artifact", "backend")
    spy <- function(name) function() {
        calls[[name]] <<- calls[[name]] + 1L
    }

    expect_error(
        auxPublicationResolveSurface(
            "unknown.surface",
            prepare = spy("prepare"),
            runner = spy("runner"),
            artifact = spy("artifact"),
            backend = spy("backend")
        ),
        class = "multischolar_publication_auxiliary_surface_excluded"
    )
    expect_identical(calls, stats::setNames(integer(4L), names(calls)))
})

test_that("auxiliary authorities are exact and remain non-promotional", {
    parameters <- .publicationAuxRecord("parameters-v1.json")
    responses <- .publicationAuxRecord("responses-v1.json")
    sources <- .publicationAuxRecord("sources-v1.json")
    splits <- .publicationAuxRecord("splits-v1.json")

    expect_silent(auxPublicationValidateParameters(parameters))
    expect_silent(auxPublicationValidateResponses(responses))
    expect_silent(auxPublicationValidateSources(sources))
    expect_silent(auxPublicationValidateSplits(splits))
    expect_length(parameters$parameters, 18L)
    expect_identical(
        vapply(
            parameters$parameters,
            `[[`,
            character(1),
            "parameter_id"
        ),
        vapply(
            parameters$consumer_registry,
            `[[`,
            character(1),
            "parameter_id"
        )
    )
    expect_length(sources$sources, 0L)
    expect_true(all(vapply(sources$decisions, function(decision) {
        identical(
            decision$cross_project_claim_status,
            "insufficient_independent_projects"
        ) && !isTRUE(decision$promotion_eligible)
    }, logical(1))))
    expect_length(splits$assignments, 20L)
    expect_length(splits$pilot_calibration_assignments, 5L)
})

.publicationAuxFixture <- function(route_id) {
    splits <- .publicationAuxRecord("splits-v1.json")
    assignment <- auxPublicationFindAssignment(
        splits,
        route_id,
        "fixture_correctness"
    )
    contract <- auxPublicationBuildContract(
        route_id,
        "fixture_correctness",
        assignment,
        NULL
    )
    payload <- auxPublicationBuildPayload(contract)
    payload_path <- tempfile(fileext = ".rds")
    auxPublicationWritePayload(payload, payload_path)
    truth <- auxPublicationBuildTruth(
        contract,
        payload,
        publicationFileDigest(payload_path)
    )
    list(
        contract = contract,
        payload = payload,
        payload_path = payload_path,
        truth = truth
    )
}

test_that("auxiliary fixture generation is byte deterministic", {
    for (route_id in names(auxPublicationRoutes())) {
        first <- .publicationAuxFixture(route_id)
        second <- .publicationAuxFixture(route_id)
        withr::defer(unlink(first$payload_path, force = TRUE))
        withr::defer(unlink(second$payload_path, force = TRUE))

        expect_identical(
            publicationFileDigest(first$payload_path),
            publicationFileDigest(second$payload_path),
            info = route_id
        )
        expect_identical(
            publicationObjectDigest(first$truth$facts),
            publicationObjectDigest(second$truth$facts),
            info = route_id
        )
    }
})

test_that("phosphosite fixture follows real exported stage graph", {
    fixture <- .publicationAuxFixture("phosphosite_stages")
    withr::defer(unlink(fixture$payload_path, force = TRUE))

    result <- suppressWarnings(auxPublicationRunPhosphosite(fixture$payload))
    expect_silent(auxPublicationValidateRunResult(
        result,
        fixture$truth,
        fixture$contract
    ))
    expect_identical(result$stages$cleaned, 20L)
    expect_identical(result$stages$probabilities, 20L)
    expect_identical(result$stages$filtered, 15L)
    expect_identical(result$stages$long, 60L)
})

test_that("multiomics fixtures execute exact local API routes", {
    expected <- list(
        mofa_weights = list(result_rows = 60L, file_count = 9L),
        metabolite_enrichment = list(result_rows = 20L, file_count = 4L),
        metabolite_pathway = list(result_rows = 10L, file_count = 1L),
        stringdb_rank = list(result_rows = 5L, file_count = 3L)
    )
    for (route_id in names(expected)) {
        fixture <- .publicationAuxFixture(route_id)
        output_root <- withr::local_tempdir()
        withr::defer(unlink(fixture$payload_path, force = TRUE))
        invisible(capture.output(
            result <- suppressWarnings(suppressMessages(
                auxPublicationRunPayload(
                    fixture$contract,
                    fixture$payload,
                    output_root
                )
            ))
        ))

        expect_silent(auxPublicationValidateRunResult(
            result,
            fixture$truth,
            fixture$contract
        ))
        result_rows <- if (identical(route_id, "mofa_weights")) {
            result$stages$plotted_feature_count
        } else {
            result$stages$result_row_count
        }
        expect_identical(result_rows, expected[[route_id]]$result_rows)
        expect_identical(
            result$stages$file_count,
            expected[[route_id]]$file_count
        )
    }
})

test_that("service response mutations fail closed", {
    responses <- .publicationAuxRecord("responses-v1.json")
    definition <- auxPublicationResponseDefinition(responses, "kegg")
    response <- auxPublicationExpandResponse(definition, 12L, 4L, 250L)

    expect_silent(auxPublicationValidateExpandedResponse(
        response,
        definition,
        12L,
        4L,
        250L
    ))

    reordered <- response[rev(seq_len(nrow(response))), , drop = FALSE]
    expect_error(
        auxPublicationValidateExpandedResponse(
            reordered,
            definition,
            12L,
            4L,
            250L
        ),
        class = "multischolar_publication_auxiliary_response_invalid"
    )

    missing <- response
    missing$mappedIDs[[1L]] <- NA_character_
    expect_error(
        auxPublicationValidateExpandedResponse(
            missing,
            definition,
            12L,
            4L,
            250L
        ),
        class = "multischolar_publication_auxiliary_response_invalid"
    )

    stale <- publicationGovernanceCopy(definition)
    stale$request_response_sha256 <- strrep("0", 64L)
    expect_error(
        auxPublicationValidateResponseDefinition(stale),
        class = "multischolar_publication_auxiliary_contract_error"
    )
})

test_that("all auxiliary negative contracts reject at admission", {
    negatives <- publicationReadJson(.publicationAuxRepoPath(
        "tests",
        "testdata",
        "omics-performance",
        "auxiliary",
        "negative",
        "contracts-v1.json"
    ))

    expect_silent(auxPublicationValidateNegatives(negatives))
    receipts <- lapply(negatives$cases, auxPublicationEvaluateNegative)
    expect_length(receipts, 22L)
    expect_true(all(vapply(receipts, function(receipt) {
        identical(receipt$network_invocations, 0L) &&
            identical(receipt$runner_invocations, 0L) &&
            identical(receipt$artifact_invocations, 0L) &&
            identical(receipt$backend_invocations, 0L) &&
            !isTRUE(receipt$promotion_authority)
    }, logical(1))))
})

test_that("auxiliary manifest fails closed into optimization handoff", {
    manifest <- .publicationAuxRecord("manifest-v1.json")

    expect_silent(auxPublicationValidateManifest(manifest))
    expect_identical(manifest$status, "optimization_required_incomplete")
    expect_length(manifest$workloads, 10L)
    expect_length(manifest$missing_workloads, 10L)
    expect_true(manifest$fixture_matrix_complete)
    expect_false(manifest$representative_matrix_complete)
    expect_false(manifest$heavy_matrix_complete)
    expect_false(manifest$stress_matrix_complete)
    expect_false(manifest$candidate_access_allowed)
    expect_false(manifest$confirmatory_benchmark_allowed)

    phosphosite <- Filter(function(route) {
        identical(route$route_id, "phosphosite_stages")
    }, manifest$route_readiness)[[1L]]
    expect_identical(
        phosphosite$representative_status,
        "historical_timeout_optimization_required"
    )
    expect_true(phosphosite$optimization_handoff_allowed)
    expect_false(phosphosite$confirmation_allowed)
})
