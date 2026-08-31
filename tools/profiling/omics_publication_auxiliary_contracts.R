.AUX_PUBLICATION_OWNER <- "OMICS-ART-068"
.AUX_PUBLICATION_VERSION <- "1.0.0"

auxPublicationAbort <- function(message, class = "contract_error") {
    publicationAbort(
        message,
        paste0("multischolar_publication_auxiliary_", class)
    )
}

auxPublicationBinding <- function(path) {
    list(path = path, sha256 = publicationFileDigest(path))
}

auxPublicationValidateBinding <- function(binding, label) {
    publicationRequireNames(binding, c("path", "sha256"), label)
    valid <- publicationScalarString(binding$path) &&
        publicationScalarString(binding$sha256) &&
        grepl("^[0-9a-f]{64}$", binding$sha256) &&
        file.exists(publicationPath(binding$path)) &&
        identical(publicationFileDigest(binding$path), binding$sha256)
    if (!valid) auxPublicationAbort(paste(label, "binding differs"))
    invisible(binding)
}

auxPublicationSurfaces <- function() {
    list(
        phosphosite = list(
            surface_id = "phosphosite.api",
            surface_type = "auxiliary_exported_api",
            primary_entry_point = "processMultisiteEvidence",
            compatibility_entry_points = list(
                "addColumnsToEvidenceTbl", "getMaxProb", "getBestPosition",
                "chooseBestPhosphositeAccession", "removePeptidesWithoutAbundances",
                "filterPeptideAndExtractProbabilities", "addPeptideStartAndEnd",
                "addPhosphositesPositionsString", "addXMerStrings",
                "filterByScoreAndGetSimilarPeptides", "allPhosphositesPivotLonger",
                "groupParalogPeptides", "allPhosphositesPivotWider",
                "uniquePhosphositesSummariseLongList",
                "uniquePhosphositesSummariseWideList"
            ),
            claim_scope = "compatibility_and_performance_only"
        ),
        multiomics = list(
            surface_id = "multiomics.api",
            surface_type = "auxiliary_exported_api",
            primary_entry_point = "runMetabolomicsPathwayEnrichment",
            compatibility_entry_points = list(
                "plotMofaWeights", "runMetabolomicsEnrichmentAnalysis",
                "runKeggEnrichment", "runReactomeEnrichment",
                "runOneStringDbRankEnrichmentMofa",
                "retrieveStringDBEnrichmentResults"
            ),
            claim_scope = "compatibility_and_performance_only"
        )
    )
}

auxPublicationValidateSurface <- function(surface) {
    publicationRequireNames(surface, c(
        "surface_id", "surface_type", "primary_entry_point",
        "compatibility_entry_points", "claim_scope"
    ), "Auxiliary surface")
    expected <- Filter(function(value) {
        identical(value$surface_id, surface$surface_id)
    }, auxPublicationSurfaces())
    valid <- length(expected) == 1L && identical(surface, expected[[1L]]) &&
        identical(surface$surface_type, "auxiliary_exported_api") &&
        identical(surface$claim_scope, "compatibility_and_performance_only")
    if (!valid) auxPublicationAbort("auxiliary surface differs")
    invisible(surface)
}

auxPublicationValidateSurfaceAuthority <- function(record) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "surfaces_id", "owner_ticket_id", "status",
        "capability_binding", "coverage_binding", "surfaces",
        "workflow_capability_creation_allowed", "artifact_state_ownership_allowed",
        "automatic_backend_decision_allowed", "unknown_policy",
        "publication_authority"
    ), "Auxiliary surface authority")
    auxPublicationValidateBinding(record$capability_binding, "Capability")
    auxPublicationValidateBinding(record$coverage_binding, "Coverage")
    lapply(record$surfaces, auxPublicationValidateSurface)
    capability <- publicationReadJson(record$capability_binding$path)
    inventory <- Filter(function(value) {
        value$surface_id %in% vapply(
            auxPublicationSurfaces(),
            `[[`,
            character(1),
            "surface_id"
        )
    }, capability$non_workflow_surfaces)
    observed <- vapply(record$surfaces, `[[`, character(1), "surface_id")
    inventory_ids <- vapply(inventory, `[[`, character(1), "surface_id")
    expected_ids <- unname(vapply(
        auxPublicationSurfaces(),
        `[[`,
        character(1),
        "surface_id"
    ))
    valid <- identical(
        record$schema,
        "multischolar.omics_publication_auxiliary_surfaces"
    ) && identical(record$schema_version, .AUX_PUBLICATION_VERSION) &&
        identical(
            record$surfaces_id,
            "multischolar.omics_publication_auxiliary_surfaces.2026-08-28.v1"
        ) &&
        identical(record$owner_ticket_id, .AUX_PUBLICATION_OWNER) &&
        identical(record$status, "frozen_pre_candidate") &&
        identical(observed, expected_ids) && !anyDuplicated(observed) &&
        setequal(observed, inventory_ids) &&
        !isTRUE(record$workflow_capability_creation_allowed) &&
        !isTRUE(record$artifact_state_ownership_allowed) &&
        !isTRUE(record$automatic_backend_decision_allowed) &&
        identical(record$unknown_policy, "reject") &&
        !isTRUE(record$publication_authority)
    if (!valid) auxPublicationAbort("auxiliary surface authority differs")
    invisible(record)
}

auxPublicationBuildSurfaceAuthority <- function() {
    list(
        schema = "multischolar.omics_publication_auxiliary_surfaces",
        schema_version = .AUX_PUBLICATION_VERSION,
        surfaces_id = paste0(
            "multischolar.omics_publication_auxiliary_surfaces.2026-08-28.v1"
        ),
        owner_ticket_id = .AUX_PUBLICATION_OWNER,
        status = "frozen_pre_candidate",
        capability_binding = auxPublicationBinding(
            "tests/testdata/omics-capabilities.json"
        ),
        coverage_binding = auxPublicationBinding(
            "tests/testdata/omics-performance/coverage-v1.json"
        ),
        surfaces = unname(auxPublicationSurfaces()),
        workflow_capability_creation_allowed = FALSE,
        artifact_state_ownership_allowed = FALSE,
        automatic_backend_decision_allowed = FALSE,
        unknown_policy = "reject",
        publication_authority = FALSE
    )
}

auxPublicationResolveSurface <- function(
    surface_id,
    prepare = function() NULL,
    runner = function() NULL,
    artifact = function() NULL,
    backend = function() NULL
) {
    ids <- vapply(auxPublicationSurfaces(), `[[`, character(1), "surface_id")
    if (!surface_id %in% ids) {
        auxPublicationAbort(
            paste("auxiliary surface is excluded:", surface_id),
            "surface_excluded"
        )
    }
    surface <- auxPublicationSurfaces()[[which(ids == surface_id)]]
    list(
        surface = surface,
        prepare = prepare,
        runner = runner,
        artifact = artifact,
        backend = backend
    )
}

auxPublicationRoutes <- function() {
    list(
        phosphosite_stages = list(
            route_id = "phosphosite_stages",
            surface_id = "phosphosite.api",
            primary_entry_point = "processMultisiteEvidence",
            measured_boundary = "staged_pipeline_after_fasta_parse",
            service_ids = list(),
            performance_work_unit = "evidence_rows_processed",
            orchestration_compatibility_binding =
                "evidence_col_to_use_existing_api_binding",
            fitting_or_inference_claim = FALSE
        ),
        mofa_weights = list(
            route_id = "mofa_weights",
            surface_id = "multiomics.api",
            primary_entry_point = "plotMofaWeights",
            measured_boundary = "weight_ranking_plot_object_and_file_products",
            service_ids = list(),
            performance_work_unit = "feature_weights_ranked",
            orchestration_compatibility_binding = "frozen_mofa_model_double",
            fitting_or_inference_claim = FALSE
        ),
        metabolite_enrichment = list(
            route_id = "metabolite_enrichment",
            surface_id = "multiomics.api",
            primary_entry_point = "runMetabolomicsEnrichmentAnalysis",
            measured_boundary = "local_rank_join_response_combine_and_write",
            service_ids = list("kegg", "reactome"),
            performance_work_unit = "feature_weights_joined",
            orchestration_compatibility_binding =
                "service_specific_immutable_local_responses",
            fitting_or_inference_claim = FALSE
        ),
        metabolite_pathway = list(
            route_id = "metabolite_pathway",
            surface_id = "multiomics.api",
            primary_entry_point = "runMetabolomicsPathwayEnrichment",
            measured_boundary = "local_assay_combine_identifier_map_and_write",
            service_ids = list("metabolomics_pathway"),
            performance_work_unit = "mapped_identifiers_processed",
            orchestration_compatibility_binding =
                "metabolomics_pathway_immutable_local_response",
            fitting_or_inference_claim = FALSE
        ),
        stringdb_rank = list(
            route_id = "stringdb_rank",
            surface_id = "multiomics.api",
            primary_entry_point = "runOneStringDbRankEnrichmentMofa",
            measured_boundary = "local_rank_response_and_file_products",
            service_ids = list("stringdb"),
            performance_work_unit = "ranked_identifiers_processed",
            orchestration_compatibility_binding =
                "stringdb_immutable_local_response",
            fitting_or_inference_claim = FALSE
        )
    )
}

auxPublicationValidateRoute <- function(route) {
    publicationRequireNames(route, c(
        "route_id", "surface_id", "primary_entry_point", "measured_boundary",
        "service_ids", "performance_work_unit",
        "orchestration_compatibility_binding", "fitting_or_inference_claim"
    ), "Auxiliary route")
    expected <- Filter(function(value) {
        identical(value$route_id, route$route_id)
    }, auxPublicationRoutes())
    valid <- length(expected) == 1L && identical(route, expected[[1L]]) &&
        !isTRUE(route$fitting_or_inference_claim)
    if (!valid) auxPublicationAbort("auxiliary route differs")
    invisible(route)
}

auxPublicationClassDimensions <- function(route_id, workload_class) {
    classes <- c(
        "fixture_correctness", "representative", "operational_heavy", "stress"
    )
    if (!workload_class %in% classes) {
        auxPublicationAbort("auxiliary workload class is unsupported")
    }
    position <- match(workload_class, classes)
    if (identical(route_id, "phosphosite_stages")) {
        rows <- c(20L, 100000L, 1000000L, 2000000L)[[position]]
        samples <- c(4L, 24L, 64L, 96L)[[position]]
        return(list(
            evidence_row_count = rows,
            feature_count = rows,
            sample_count = samples,
            measured_sample_count = samples,
            layer_count = 1L,
            response_row_count = 0L,
            primary_work_units = as.numeric(rows),
            sample_scaling_measured = TRUE
        ))
    }
    features <- c(1000L, 50000L, 100000L, 250000L)[[position]]
    samples <- c(20L, 100L, 200L, 400L)[[position]]
    response_rows <- if (route_id %in% c(
        "metabolite_enrichment",
        "metabolite_pathway",
        "stringdb_rank"
    )) {
        c(5L, 10000L, 25000L, 50000L)[[position]]
    } else {
        0L
    }
    list(
        evidence_row_count = 0L,
        feature_count = features,
        sample_count = samples,
        measured_sample_count = if (route_id %in% c(
            "metabolite_enrichment",
            "metabolite_pathway"
        )) 1L else 0L,
        layer_count = 3L,
        response_row_count = response_rows,
        primary_work_units = as.numeric(features),
        sample_scaling_measured = FALSE
    )
}

auxPublicationValidateDimensions <- function(dimensions, route_id) {
    publicationRequireNames(dimensions, c(
        "evidence_row_count", "feature_count", "sample_count",
        "measured_sample_count", "layer_count", "response_row_count",
        "primary_work_units", "sample_scaling_measured"
    ), "Auxiliary dimensions")
    numeric_fields <- unlist(dimensions[c(
        "evidence_row_count", "feature_count", "sample_count",
        "measured_sample_count", "layer_count", "response_row_count",
        "primary_work_units"
    )], use.names = FALSE)
    valid <- all(is.finite(numeric_fields)) && all(numeric_fields >= 0) &&
        dimensions$feature_count > 0 && dimensions$sample_count > 0 &&
        dimensions$primary_work_units > 0 &&
        identical(
            isTRUE(dimensions$sample_scaling_measured),
            identical(route_id, "phosphosite_stages")
        )
    if (!valid) auxPublicationAbort("auxiliary dimensions differ")
    invisible(dimensions)
}

auxPublicationWorkloadId <- function(route_id, workload_class) {
    paste("auxiliary", route_id, workload_class, "v1", sep = ".")
}

auxPublicationContractSourceFiles <- function() {
    list(
        "tools/profiling/omics_publication_auxiliary_contracts.R",
        "tools/profiling/omics_publication_auxiliary_model.R",
        "tools/profiling/omics_publication_auxiliary_truth.R",
        "tools/profiling/omics_publication_workload_auxiliary.R"
    )
}

auxPublicationBuildContract <- function(
    route_id,
    workload_class,
    assignment,
    seed,
    expected_payload_sha256 = strrep("0", 64L),
    expected_truth_sha256 = strrep("0", 64L)
) {
    route <- auxPublicationRoutes()[[route_id]]
    if (is.null(route)) auxPublicationAbort("auxiliary route is unsupported")
    dimensions <- auxPublicationClassDimensions(route_id, workload_class)
    list(
        schema = "multischolar.omics_publication_auxiliary_workload",
        schema_version = .AUX_PUBLICATION_VERSION,
        contract_id = paste0(
            auxPublicationWorkloadId(route_id, workload_class),
            ".contract"
        ),
        owner_ticket_id = .AUX_PUBLICATION_OWNER,
        status = "frozen_pre_candidate",
        workload_id = auxPublicationWorkloadId(route_id, workload_class),
        workload_class = workload_class,
        route = route,
        dimensions = dimensions,
        parameter_authority = auxPublicationBinding(paste0(
            "tests/testdata/omics-performance/auxiliary/parameters-v1.json"
        )),
        source_authority = auxPublicationBinding(paste0(
            "tests/testdata/omics-performance/auxiliary/sources-v1.json"
        )),
        response_authority = auxPublicationBinding(paste0(
            "tests/testdata/omics-performance/auxiliary/responses-v1.json"
        )),
        split_authority = auxPublicationBinding(paste0(
            "tests/testdata/omics-performance/auxiliary/splits-v1.json"
        )),
        surface_authority = auxPublicationBinding(paste0(
            "tests/testdata/omics-performance/auxiliary/surfaces-v1.json"
        )),
        assignment = assignment,
        generator = list(
            mode = "repository_generated",
            source_bindings = lapply(
                auxPublicationContractSourceFiles(),
                auxPublicationBinding
            ),
            payload_filename = "payload.rds",
            truth_filename = "truth.json",
            generated_outside_measurement = TRUE
        ),
        rng = list(
            kind = "L'Ecuyer-CMRG",
            normal_kind = "Inversion",
            sample_kind = "Rejection",
            seed_family_id = assignment$seed_family_id,
            seed = seed
        ),
        execution = list(
            process_isolation = "fresh_vanilla_R",
            candidate_loaded = FALSE,
            historical_pilot_required = identical(
                workload_class,
                "operational_heavy"
            ),
            swap_allowed = FALSE,
            timeout_seconds = 900L,
            network_allowed = FALSE
        ),
        truth_contract = list(
            schema_id = "auxiliary.truth.v1",
            oracle_method = "independent_generated_truth",
            exact_order_required = TRUE,
            service_behaviour_authority = FALSE,
            biological_authority = FALSE,
            mofa_fitting_authority = FALSE
        ),
        claim_scope = list(
            evidence_class = "public_generated",
            compatibility_authority = TRUE,
            local_performance_authority = TRUE,
            cross_project_authority = FALSE,
            workflow_authority = FALSE,
            artifact_authority = FALSE,
            backend_policy_authority = FALSE,
            biological_authority = FALSE,
            promotion_authority = FALSE
        ),
        expected_digests = list(
            payload_sha256 = expected_payload_sha256,
            truth_sha256 = expected_truth_sha256
        ),
        publication_authority = FALSE
    )
}

auxPublicationValidateContract <- function(contract) {
    publicationRequireNames(contract, c(
        "schema", "schema_version", "contract_id", "owner_ticket_id", "status",
        "workload_id", "workload_class", "route", "dimensions",
        "parameter_authority", "source_authority", "response_authority",
        "split_authority", "surface_authority", "assignment", "generator",
        "rng", "execution", "truth_contract", "claim_scope",
        "expected_digests", "publication_authority"
    ), "Auxiliary workload contract")
    auxPublicationValidateRoute(contract$route)
    auxPublicationValidateDimensions(contract$dimensions, contract$route$route_id)
    lapply(c(
        "parameter_authority", "source_authority", "response_authority",
        "split_authority", "surface_authority"
    ), function(name) {
        auxPublicationValidateBinding(contract[[name]], name)
    })
    publicationRequireNames(contract$expected_digests, c(
        "payload_sha256", "truth_sha256"
    ), "Auxiliary expected digests")
    publicationRequireNames(contract$generator, c(
        "mode", "source_bindings", "payload_filename", "truth_filename",
        "generated_outside_measurement"
    ), "Auxiliary generator")
    publicationRequireNames(contract$rng, c(
        "kind", "normal_kind", "sample_kind", "seed_family_id", "seed"
    ), "Auxiliary RNG")
    publicationRequireNames(contract$execution, c(
        "process_isolation", "candidate_loaded", "historical_pilot_required",
        "swap_allowed", "timeout_seconds", "network_allowed"
    ), "Auxiliary execution")
    publicationRequireNames(contract$truth_contract, c(
        "schema_id", "oracle_method", "exact_order_required",
        "service_behaviour_authority", "biological_authority",
        "mofa_fitting_authority"
    ), "Auxiliary truth contract")
    publicationRequireNames(contract$claim_scope, c(
        "evidence_class", "compatibility_authority",
        "local_performance_authority", "cross_project_authority",
        "workflow_authority", "artifact_authority",
        "backend_policy_authority", "biological_authority",
        "promotion_authority"
    ), "Auxiliary claim scope")
    splits <- publicationReadJson(contract$split_authority$path)
    if (exists("auxPublicationValidateSplits", mode = "function")) {
        auxPublicationValidateSplits(splits)
    }
    expected_assignment <- Filter(function(assignment) {
        identical(assignment$route_id, contract$route$route_id) &&
            identical(assignment$workload_class, contract$workload_class) &&
            !identical(assignment$split_role, "pilot")
    }, splits$assignments)
    expected_source_paths <- unlist(
        auxPublicationContractSourceFiles(),
        use.names = FALSE
    )
    source_paths <- vapply(
        contract$generator$source_bindings,
        `[[`,
        character(1),
        "path"
    )
    source_bindings_valid <- identical(source_paths, expected_source_paths) &&
        !anyDuplicated(source_paths) && all(vapply(
            contract$generator$source_bindings,
            function(binding) {
                tryCatch({
                    auxPublicationValidateBinding(binding, "Generator source")
                    TRUE
                }, error = function(error) FALSE)
            },
            logical(1)
        ))
    digests <- unlist(contract$expected_digests, use.names = FALSE)
    valid <- identical(
        contract$schema,
        "multischolar.omics_publication_auxiliary_workload"
    ) && identical(contract$schema_version, .AUX_PUBLICATION_VERSION) &&
        identical(contract$owner_ticket_id, .AUX_PUBLICATION_OWNER) &&
        identical(contract$status, "frozen_pre_candidate") &&
        identical(
            contract$contract_id,
            paste0(contract$workload_id, ".contract")
        ) &&
        identical(
            contract$workload_id,
            auxPublicationWorkloadId(
                contract$route$route_id,
                contract$workload_class
            )
        ) &&
        identical(
            publicationObjectDigest(contract$dimensions),
            publicationObjectDigest(auxPublicationClassDimensions(
                contract$route$route_id,
                contract$workload_class
            ))
        ) &&
        length(expected_assignment) == 1L && identical(
            publicationObjectDigest(contract$assignment),
            publicationObjectDigest(expected_assignment[[1L]])
        ) &&
        source_bindings_valid &&
        identical(contract$generator$mode, "repository_generated") &&
        identical(contract$generator$payload_filename, "payload.rds") &&
        identical(contract$generator$truth_filename, "truth.json") &&
        isTRUE(contract$generator$generated_outside_measurement) &&
        identical(contract$rng$kind, "L'Ecuyer-CMRG") &&
        identical(contract$rng$normal_kind, "Inversion") &&
        identical(contract$rng$sample_kind, "Rejection") &&
        identical(
            contract$rng$seed_family_id,
            contract$assignment$seed_family_id
        ) && identical(contract$rng$seed, contract$assignment$seed) &&
        all(grepl("^[0-9a-f]{64}$", digests)) &&
        identical(contract$execution$process_isolation, "fresh_vanilla_R") &&
        identical(
            contract$execution$historical_pilot_required,
            identical(contract$workload_class, "operational_heavy")
        ) &&
        identical(contract$execution$timeout_seconds, 900L) &&
        !isTRUE(contract$execution$candidate_loaded) &&
        !isTRUE(contract$execution$swap_allowed) &&
        !isTRUE(contract$execution$network_allowed) &&
        identical(contract$truth_contract$schema_id, "auxiliary.truth.v1") &&
        identical(
            contract$truth_contract$oracle_method,
            "independent_generated_truth"
        ) && isTRUE(contract$truth_contract$exact_order_required) &&
        !isTRUE(contract$truth_contract$service_behaviour_authority) &&
        !isTRUE(contract$truth_contract$biological_authority) &&
        !isTRUE(contract$truth_contract$mofa_fitting_authority) &&
        identical(contract$claim_scope$evidence_class, "public_generated") &&
        isTRUE(contract$claim_scope$compatibility_authority) &&
        isTRUE(contract$claim_scope$local_performance_authority) &&
        !isTRUE(contract$claim_scope$cross_project_authority) &&
        !isTRUE(contract$claim_scope$workflow_authority) &&
        !isTRUE(contract$claim_scope$artifact_authority) &&
        !isTRUE(contract$claim_scope$backend_policy_authority) &&
        !isTRUE(contract$claim_scope$biological_authority) &&
        !isTRUE(contract$claim_scope$promotion_authority) &&
        !isTRUE(contract$publication_authority)
    if (!valid) auxPublicationAbort("auxiliary workload contract differs")
    invisible(contract)
}
