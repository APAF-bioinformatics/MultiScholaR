metabPublicationAuthorityPaths <- function() {
    root <- "tests/testdata/omics-performance/metabolomics"
    list(
        parameters = file.path(root, "parameters-v1.json"),
        sources = file.path(root, "sources-v1.json"),
        splits = file.path(root, "splits-v1.json"),
        mappings = file.path(root, "mapping-authorities-v1.json"),
        boundaries = file.path(root, "support-boundaries-v1.json")
    )
}

metabPublicationAuthorityBinding <- function(path) {
    list(path = path, sha256 = publicationFileDigest(path))
}

metabPublicationAssayCounts <- function(profile_id, aggregate_count) {
    assays <- metabPublicationFixtureAssays(profile_id)
    base <- aggregate_count %/% length(assays)
    remainder <- aggregate_count %% length(assays)
    counts <- rep(base, length(assays))
    if (remainder) counts[seq_len(remainder)] <- counts[seq_len(remainder)] + 1L
    stats::setNames(as.list(as.integer(counts)), assays)
}

metabPublicationGeneratedDimensions <- function(
    route,
    profile_id,
    workload_class
) {
    scale <- metabPublicationExpectedScale(profile_id, workload_class)
    aggregate <- scale$minimum_features
    samples <- scale$sample_count
    counts <- metabPublicationAssayCounts(profile_id, aggregate)
    list(
        aggregate_feature_count = as.integer(aggregate),
        assay_feature_counts = counts,
        sample_count = as.integer(samples),
        assay_count = as.integer(length(counts)),
        quantity_count = as.numeric(aggregate * samples),
        payload_member_count = if (identical(route, "msdial") &&
            identical(profile_id, "mixed")) 3L else 1L
    )
}

metabPublicationSplitAssignmentFor <- function(
    splits,
    route,
    profile_id,
    workload_class
) {
    workload_id <- metabPublicationWorkloadId(route, profile_id, workload_class)
    found <- Filter(
        function(value) identical(value$assignment_id, workload_id),
        splits$assignments
    )
    if (length(found) != 1L) {
        metabPublicationAbort("metabolomics workload split is absent")
    }
    found[[1L]]
}

metabPublicationGeneratorSources <- function() {
    paths <- c(
        "tools/profiling/omics_publication_protocol.R",
        "tools/profiling/omics_publication_metabolomics_contracts.R",
        "tools/profiling/omics_publication_metabolomics_model.R",
        "tools/profiling/omics_publication_metabolomics_serializers.R",
        "tools/profiling/omics_publication_metabolomics_truth.R",
        "tools/profiling/omics_publication_metabolomics_sources.R",
        "tools/profiling/omics_publication_metabolomics_governance.R",
        "tools/profiling/omics_publication_metabolomics_fixtures.R",
        "tools/profiling/omics_publication_workload_metabolomics.R"
    )
    lapply(as.list(paths), metabPublicationAuthorityBinding)
}

metabPublicationRng <- function(assignment, workload_class) {
    fixture <- identical(workload_class, "fixture_correctness")
    seed <- assignment$seed
    streams <- if (fixture) {
        list()
    } else {
        as.list(stats::setNames(
            as.integer(seed + seq_len(7L) * 100L),
            c(
                "hierarchy", "feature_offsets", "batch", "residual",
                "mcar", "mar", "mnar"
            )
        ))
    }
    list(
        kind = "L'Ecuyer-CMRG",
        normal_kind = "Inversion",
        sample_kind = "Rejection",
        seed_family_id = assignment$seed_family_id,
        seed = seed,
        streams = streams
    )
}

metabPublicationFixtureGenerator <- function(route, profile_id, fixture_root) {
    root <- file.path(fixture_root, paste(route, profile_id, sep = "-"))
    members <- unlist(
        metabPublicationOutputMembers(route, profile_id),
        use.names = FALSE
    )
    payload_paths <- file.path(root, "payload", members)
    source_paths <- c(
        metabPublicationFixtureSourcePaths(route, profile_id),
        metabPublicationFixtureReviewPaths(profile_id),
        "tools/profiling/omics_publication_metabolomics_fixtures.R"
    )
    list(
        mode = "fixture_replay",
        source_bindings = lapply(
            as.list(unique(source_paths)),
            metabPublicationAuthorityBinding
        ),
        chunk_rows = NULL,
        output_members = as.list(members),
        truth_filename = "truth.json",
        fixture_payloads = lapply(
            as.list(payload_paths),
            metabPublicationAuthorityBinding
        ),
        fixture_truth = metabPublicationAuthorityBinding(
            file.path(root, "truth.json")
        )
    )
}

metabPublicationGeneratedGenerator <- function(route, profile_id) {
    list(
        mode = "generated",
        source_bindings = metabPublicationGeneratorSources(),
        chunk_rows = 5000L,
        output_members = metabPublicationOutputMembers(route, profile_id),
        truth_filename = "truth.json",
        fixture_payloads = list(),
        fixture_truth = NULL
    )
}

metabPublicationTruthContract <- function(workload_class) {
    fixture <- identical(workload_class, "fixture_correctness")
    list(
        schema_id = "metabolomics.truth.v1",
        oracle_method = if (fixture) {
            "direct_raw_table_arithmetic_and_reviewed_e2e_authority"
        } else {
            "independent_online_generated_truth"
        },
        validated_domains = as.list(c(
            "payload", "mapping", "assay_identity", "dimensions", "quantities",
            "missingness", "effects", "design"
        )),
        independent_of_package_reader = TRUE,
        support_widening_allowed = FALSE
    )
}

metabPublicationExecution <- function(dimensions, workload_class) {
    list(
        preparation_processes = 2L,
        process_isolation = "fresh_vanilla_R",
        candidate_loaded = FALSE,
        historical_pilot_required = identical(
            workload_class,
            "operational_heavy"
        ),
        swap_allowed = FALSE,
        timeout_seconds = if (workload_class %in% c(
            "operational_heavy",
            "stress"
        )) 1800 else 900,
        expected_peak_bytes = as.numeric(dimensions$quantity_count * 16)
    )
}

metabPublicationPrivacy <- function(workload_class) {
    list(
        classification = if (identical(
            workload_class,
            "fixture_correctness"
        )) "public_fixture" else "public_generated",
        private_source_opened = FALSE,
        payload_tracked = FALSE,
        direct_identifiers_present = FALSE,
        cross_omic_values_reused = FALSE
    )
}

metabPublicationClaimScope <- function(route, workload_class) {
    fixture <- identical(workload_class, "fixture_correctness")
    performance <- workload_class %in% c(
        "representative",
        "operational_heavy"
    )
    list(
        evidence_class = switch(
            workload_class,
            fixture_correctness = "fixture_correctness",
            stress = "stress_characterization",
            "public_generated"
        ),
        verified_stages = metabPublicationCapabilities()[[route]]$verified_stages,
        scientific_authority = fixture,
        performance_authority = performance,
        cross_project_authority = FALSE,
        promotion_authority = FALSE,
        limitations = list(
            "Project-specific evidence without independent real-project sources.",
            "Generated evidence cannot widen the declared route support boundary."
        )
    )
}

metabPublicationBuildContract <- function(
    route,
    profile_id,
    workload_class,
    fixture_root = "tests/testdata/omics-performance/metabolomics/fixtures"
) {
    paths <- metabPublicationAuthorityPaths()
    splits <- publicationReadJson(paths$splits)
    assignment <- metabPublicationSplitAssignmentFor(
        splits,
        route,
        profile_id,
        workload_class
    )
    fixture <- identical(workload_class, "fixture_correctness")
    fixture_truth <- if (fixture) {
        publicationReadJson(file.path(
            fixture_root,
            paste(route, profile_id, sep = "-"),
            "truth.json"
        ))
    } else {
        NULL
    }
    dimensions <- if (fixture) {
        fixture_truth$dimensions
    } else {
        metabPublicationGeneratedDimensions(route, profile_id, workload_class)
    }
    member_count <- dimensions$payload_member_count
    contract <- list(
        schema = "multischolar.omics_publication_metabolomics_workload",
        schema_version = .METAB_PUBLICATION_VERSION,
        contract_id = paste(
            "metabolomics", route, profile_id, workload_class, "v1", sep = "."
        ),
        owner_ticket_id = .METAB_PUBLICATION_OWNER,
        status = "frozen",
        workload_id = assignment$assignment_id,
        workload_class = workload_class,
        capability = metabPublicationExpectedCapability(route),
        assay_profile = list(
            profile_id = profile_id,
            assays = metabPublicationAssayProfiles()[[profile_id]]$assays,
            payload_mode = if (member_count == 3L) "bundle" else "single",
            member_count = as.integer(member_count)
        ),
        dimensions = dimensions,
        model_profile_id = assignment$generator_parameter_profile_id,
        parameter_authority = metabPublicationAuthorityBinding(paths$parameters),
        source_authority = metabPublicationAuthorityBinding(paths$sources),
        split_authority = metabPublicationAuthorityBinding(paths$splits),
        mapping_authority = metabPublicationAuthorityBinding(paths$mappings),
        contract_mapping_id = paste(
            "metabolomics.mapping",
            route,
            profile_id,
            sep = "."
        ),
        support_boundary = metabPublicationAuthorityBinding(paths$boundaries),
        generator = if (fixture) {
            metabPublicationFixtureGenerator(route, profile_id, fixture_root)
        } else {
            metabPublicationGeneratedGenerator(route, profile_id)
        },
        rng = metabPublicationRng(assignment, workload_class),
        truth_contract = metabPublicationTruthContract(workload_class),
        execution = metabPublicationExecution(dimensions, workload_class),
        privacy = metabPublicationPrivacy(workload_class),
        claim_scope = metabPublicationClaimScope(route, workload_class),
        expected_digests = list(
            payload_set_sha256 = if (fixture) {
                fixture_truth$payload$payload_set_sha256
            } else {
                strrep("0", 64L)
            },
            truth_sha256 = if (fixture) {
                publicationFileDigest(file.path(
                    fixture_root,
                    paste(route, profile_id, sep = "-"),
                    "truth.json"
                ))
            } else {
                strrep("0", 64L)
            }
        ),
        publication_authority = FALSE
    )
    metabPublicationValidateWorkload(contract)
    contract
}

metabPublicationContractFilename <- function(
    route,
    profile_id,
    workload_class
) {
    paste(
        gsub("_", "-", workload_class),
        route,
        gsub("_", "-", profile_id),
        "v1.json",
        sep = "-"
    )
}

metabPublicationFreezeFixtureContracts <- function(output_root) {
    dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
    paths <- character()
    for (route in names(metabPublicationCapabilities())) {
        for (profile_id in names(metabPublicationAssayProfiles())) {
            contract <- metabPublicationBuildContract(
                route,
                profile_id,
                "fixture_correctness"
            )
            path <- file.path(output_root, metabPublicationContractFilename(
                route,
                profile_id,
                "fixture_correctness"
            ))
            publicationWriteJson(contract, path)
            paths <- c(paths, path)
        }
    }
    paths
}
