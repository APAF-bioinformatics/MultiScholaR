auxPublicationManifestWorkloadFiles <- function() {
    root <- "tests/testdata/omics-performance/auxiliary/workloads"
    paths <- sort(list.files(
        publicationPath(root),
        pattern = "[.]json$",
        full.names = FALSE
    ), method = "radix")
    as.list(file.path(root, paths))
}

auxPublicationManifestEntries <- function() {
    lapply(auxPublicationManifestWorkloadFiles(), function(path) {
        contract <- auxPublicationReadContract(path)
        truth_path <- file.path(
            "tests/testdata/omics-performance/auxiliary/truth",
            paste0(contract$workload_id, ".json")
        )
        truth <- publicationReadJson(truth_path)
        list(
            workload_id = contract$workload_id,
            route_id = contract$route$route_id,
            surface_id = contract$route$surface_id,
            workload_class = contract$workload_class,
            contract = list(
                path = path,
                sha256 = publicationFileDigest(path)
            ),
            truth = list(
                path = truth_path,
                sha256 = publicationFileDigest(truth_path)
            ),
            payload_sha256 = contract$expected_digests$payload_sha256,
            dimensions = contract$dimensions,
            exact_truth_validated = identical(
                truth$contract_basis_sha256,
                auxPublicationContractBasis(contract)
            ),
            candidate_loaded = FALSE,
            promotion_authority = FALSE
        )
    })
}

auxPublicationRouteReadiness <- function() {
    list(
        list(
            route_id = "phosphosite_stages",
            fixture_status = "exact_comparator_parity",
            representative_status = "candidate_exact_historical_timeout_resolved",
            operational_heavy_status = "frozen_host_capacity_blocked",
            stress_status = "generated_characterization_only",
            optimization_handoff_allowed = FALSE,
            confirmation_allowed = FALSE,
            promotion_authority = FALSE
        ),
        list(
            route_id = "mofa_weights",
            fixture_status = "exact_comparator_parity",
            representative_status = "exact_comparator_parity",
            operational_heavy_status = "candidate_exact_truth",
            stress_status = "generated_characterization_only",
            optimization_handoff_allowed = FALSE,
            confirmation_allowed = FALSE,
            promotion_authority = FALSE
        ),
        list(
            route_id = "metabolite_enrichment",
            fixture_status = "exact_comparator_parity",
            representative_status = "exact_comparator_parity",
            operational_heavy_status = "candidate_exact_truth",
            stress_status = "generated_characterization_only",
            optimization_handoff_allowed = FALSE,
            confirmation_allowed = FALSE,
            promotion_authority = FALSE
        ),
        list(
            route_id = "metabolite_pathway",
            fixture_status = "exact_comparator_parity",
            representative_status = "exact_candidate_historical_parity",
            operational_heavy_status = "candidate_exact_truth",
            stress_status = "generated_characterization_only",
            optimization_handoff_allowed = FALSE,
            confirmation_allowed = FALSE,
            promotion_authority = FALSE
        ),
        list(
            route_id = "stringdb_rank",
            fixture_status = "exact_comparator_parity",
            representative_status = "exact_comparator_parity",
            operational_heavy_status = "candidate_exact_truth",
            stress_status = "generated_characterization_only",
            optimization_handoff_allowed = FALSE,
            confirmation_allowed = FALSE,
            promotion_authority = FALSE
        )
    )
}

auxPublicationMissingWorkloads <- function(entries) {
    expected <- unlist(lapply(names(auxPublicationRoutes()), function(route_id) {
        vapply(c("operational_heavy", "stress"), function(workload_class) {
            auxPublicationWorkloadId(route_id, workload_class)
        }, character(1))
    }), use.names = FALSE)
    observed <- vapply(entries, `[[`, character(1), "workload_id")
    as.list(setdiff(expected, observed))
}

auxPublicationBuildManifest <- function() {
    entries <- auxPublicationManifestEntries()
    list(
        schema = "multischolar.omics_publication_auxiliary_manifest",
        schema_version = .AUX_PUBLICATION_VERSION,
        manifest_id = paste0(
            "multischolar.omics_publication_auxiliary_manifest.2026-08-30.v1"
        ),
        owner_ticket_id = .AUX_PUBLICATION_OWNER,
        status = "workload_matrix_complete_candidate_freeze_blocked",
        surfaces = auxPublicationBinding(paste0(
            "tests/testdata/omics-performance/auxiliary/surfaces-v1.json"
        )),
        parameters = auxPublicationBinding(paste0(
            "tests/testdata/omics-performance/auxiliary/parameters-v1.json"
        )),
        sources = auxPublicationBinding(paste0(
            "tests/testdata/omics-performance/auxiliary/sources-v1.json"
        )),
        responses = auxPublicationBinding(paste0(
            "tests/testdata/omics-performance/auxiliary/responses-v1.json"
        )),
        splits = auxPublicationBinding(paste0(
            "tests/testdata/omics-performance/auxiliary/splits-v1.json"
        )),
        exclusions = auxPublicationBinding(paste0(
            "tests/testdata/omics-performance/auxiliary/exclusions-v1.json"
        )),
        negatives = auxPublicationBinding(paste0(
            "tests/testdata/omics-performance/auxiliary/negative/",
            "contracts-v1.json"
        )),
        workloads = entries,
        route_readiness = auxPublicationRouteReadiness(),
        missing_workloads = auxPublicationMissingWorkloads(entries),
        source_minimum_satisfied = FALSE,
        fixture_matrix_complete = TRUE,
        representative_matrix_complete = TRUE,
        heavy_matrix_complete = TRUE,
        stress_matrix_complete = TRUE,
        candidate_access_allowed = FALSE,
        confirmatory_benchmark_allowed = FALSE,
        publication_authority = FALSE
    )
}

auxPublicationValidateManifest <- function(record) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "manifest_id", "owner_ticket_id", "status",
        "surfaces", "parameters", "sources", "responses", "splits",
        "exclusions", "negatives", "workloads", "route_readiness",
        "missing_workloads", "source_minimum_satisfied",
        "fixture_matrix_complete", "representative_matrix_complete",
        "heavy_matrix_complete", "stress_matrix_complete",
        "candidate_access_allowed", "confirmatory_benchmark_allowed",
        "publication_authority"
    ), "Auxiliary manifest")
    lapply(c(
        "surfaces", "parameters", "sources", "responses", "splits",
        "exclusions", "negatives"
    ), function(name) {
        auxPublicationValidateBinding(record[[name]], name)
    })
    expected_entries <- auxPublicationManifestEntries()
    ids <- vapply(record$workloads, `[[`, character(1), "workload_id")
    readiness_ids <- vapply(
        record$route_readiness,
        `[[`,
        character(1),
        "route_id"
    )
    entries_valid <- all(vapply(record$workloads, function(entry) {
        publicationRequireNames(entry, c(
            "workload_id", "route_id", "surface_id", "workload_class",
            "contract", "truth", "payload_sha256", "dimensions",
            "exact_truth_validated", "candidate_loaded", "promotion_authority"
        ), "Auxiliary manifest workload")
        auxPublicationValidateBinding(entry$contract, "Manifest contract")
        auxPublicationValidateBinding(entry$truth, "Manifest truth")
        isTRUE(entry$exact_truth_validated) &&
            !isTRUE(entry$candidate_loaded) &&
            !isTRUE(entry$promotion_authority)
    }, logical(1)))
    valid <- identical(
        record$schema,
        "multischolar.omics_publication_auxiliary_manifest"
    ) && identical(record$schema_version, .AUX_PUBLICATION_VERSION) &&
        identical(record$owner_ticket_id, .AUX_PUBLICATION_OWNER) &&
        identical(
            record$status,
            "workload_matrix_complete_candidate_freeze_blocked"
        ) &&
        identical(
            publicationObjectDigest(record$workloads),
            publicationObjectDigest(expected_entries)
        ) && length(ids) == 20L && !anyDuplicated(ids) && entries_valid &&
        identical(readiness_ids, names(auxPublicationRoutes())) &&
        identical(
            record$missing_workloads,
            auxPublicationMissingWorkloads(record$workloads)
        ) && !length(record$missing_workloads) &&
        !isTRUE(record$source_minimum_satisfied) &&
        isTRUE(record$fixture_matrix_complete) &&
        isTRUE(record$representative_matrix_complete) &&
        isTRUE(record$heavy_matrix_complete) &&
        isTRUE(record$stress_matrix_complete) &&
        !isTRUE(record$candidate_access_allowed) &&
        !isTRUE(record$confirmatory_benchmark_allowed) &&
        !isTRUE(record$publication_authority)
    if (!valid) auxPublicationAbort("auxiliary manifest differs")
    invisible(record)
}

auxPublicationBuildOptimizationHandoff <- function(manifest_path) {
    manifest <- publicationReadJson(manifest_path)
    auxPublicationValidateManifest(manifest)
    list(
        schema = "multischolar.omics_publication_auxiliary_handoff",
        schema_version = .AUX_PUBLICATION_VERSION,
        owner_ticket_id = .AUX_PUBLICATION_OWNER,
        status = "workload_matrix_complete_candidate_freeze_blocked",
        manifest = list(
            path = manifest_path,
            sha256 = publicationFileDigest(manifest_path)
        ),
        consumer_ticket_id = "OMICS-ART-077",
        confirmatory_consumer_ticket_id = "OMICS-ART-082",
        optimization_inputs_allowed = FALSE,
        confirmatory_inputs_allowed = FALSE,
        phosphosite_required_outcome = paste0(
            "complete_100k_representative_within_900_seconds_passed"
        ),
        workload_floor_reduction_allowed = FALSE,
        source_minimum_satisfied = FALSE,
        candidate_access_allowed = FALSE,
        promotion_authority = FALSE,
        publication_authority = FALSE
    )
}
