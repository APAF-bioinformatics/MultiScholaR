lipidPublicationAuthorityPaths <- function() {
    root <- "tests/testdata/omics-performance/lipidomics"
    list(
        parameters = file.path(root, "parameters-v1.json"),
        sources = file.path(root, "sources-v1.json"),
        splits = file.path(root, "splits-v1.json"),
        mappings = file.path(root, "mapping-authorities-v1.json"),
        boundaries = file.path(root, "support-boundaries-v1.json")
    )
}

lipidPublicationAuthorityBinding <- function(path) {
    list(path = path, sha256 = publicationFileDigest(path))
}

lipidPublicationAssayCounts <- function(profile_id, aggregate_count) {
    assays <- lipidPublicationFixtureAssays(profile_id)
    base <- aggregate_count %/% length(assays)
    remainder <- aggregate_count %% length(assays)
    counts <- rep(base, length(assays))
    if (remainder) counts[seq_len(remainder)] <- counts[seq_len(remainder)] + 1L
    stats::setNames(as.list(as.integer(counts)), assays)
}

lipidPublicationGeneratedDimensions <- function(
    route,
    profile_id,
    workload_class
) {
    scale <- lipidPublicationExpectedScale(profile_id, workload_class)
    aggregate <- scale$minimum_features
    samples <- scale$sample_count
    counts <- lipidPublicationAssayCounts(profile_id, aggregate)
    list(
        aggregate_feature_count = as.integer(aggregate),
        assay_feature_counts = counts,
        sample_count = as.integer(samples),
        assay_count = as.integer(length(counts)),
        quantity_count = as.numeric(aggregate * samples),
        payload_member_count = if (startsWith(profile_id, "mixed_")) 2L else 1L
    )
}

lipidPublicationSplitAssignmentFor <- function(
    splits,
    route,
    profile_id,
    workload_class
) {
    workload_id <- lipidPublicationWorkloadId(route, profile_id, workload_class)
    found <- Filter(
        function(value) identical(value$assignment_id, workload_id),
        splits$assignments
    )
    if (length(found) != 1L) {
        lipidPublicationAbort("lipidomics workload split is absent")
    }
    found[[1L]]
}

lipidPublicationGeneratorSources <- function() {
    paths <- c(
        "tools/profiling/omics_publication_protocol.R",
        "tools/profiling/omics_publication_lipidomics_contracts.R",
        "tools/profiling/omics_publication_lipidomics_model.R",
        "tools/profiling/omics_publication_lipidomics_serializers.R",
        "tools/profiling/omics_publication_lipidomics_truth.R",
        "tools/profiling/omics_publication_lipidomics_authorities.R",
        "tools/profiling/omics_publication_lipidomics_sources.R",
        "tools/profiling/omics_publication_lipidomics_governance.R",
        "tools/profiling/omics_publication_lipidomics_fixtures.R",
        "tools/profiling/omics_publication_workload_lipidomics.R"
    )
    lapply(as.list(paths), lipidPublicationAuthorityBinding)
}

lipidPublicationRng <- function(assignment, workload_class) {
    fixture <- identical(workload_class, "fixture_correctness")
    seed <- assignment$seed
    streams <- if (fixture) {
        list()
    } else {
        as.list(stats::setNames(
            as.integer(seed + seq_len(8L) * 100L),
            c(
                "hierarchy", "feature_offsets", "class_residual", "batch",
                "residual", "mcar", "mar", "mnar"
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

lipidPublicationFixtureGenerator <- function(route, profile_id, fixture_root) {
    root <- file.path(fixture_root, paste(route, profile_id, sep = "-"))
    members <- unlist(
        lipidPublicationOutputMembers(route, profile_id),
        use.names = FALSE
    )
    payload_paths <- file.path(root, "payload", members)
    truth_record <- publicationReadJson(file.path(root, "truth.json"))
    source_paths <- vapply(
        truth_record$source_authority_bindings,
        `[[`,
        character(1),
        "path"
    )
    list(
        mode = "fixture_replay",
        source_bindings = lapply(
            as.list(unique(source_paths)),
            lipidPublicationAuthorityBinding
        ),
        chunk_rows = NULL,
        output_members = as.list(members),
        truth_filename = "truth.json",
        fixture_payloads = lapply(
            as.list(payload_paths),
            lipidPublicationAuthorityBinding
        ),
        fixture_truth = lipidPublicationAuthorityBinding(
            file.path(root, "truth.json")
        )
    )
}

lipidPublicationGeneratedGenerator <- function(route, profile_id) {
    list(
        mode = "generated",
        source_bindings = lipidPublicationGeneratorSources(),
        chunk_rows = 5000L,
        output_members = lipidPublicationOutputMembers(route, profile_id),
        truth_filename = "truth.json",
        fixture_payloads = list(),
        fixture_truth = NULL
    )
}

lipidPublicationTruthContract <- function(route, profile_id, workload_class) {
    fixture <- identical(workload_class, "fixture_correctness")
    list(
        schema_id = "lipidomics.truth.v1",
        oracle_method = if (fixture) {
            lipidPublicationFixtureOracleMethod(route, profile_id)
        } else {
            "independent_online_generated_truth"
        },
        validated_domains = as.list(c(
            "payload", "mapping", "assay_identity", "dimensions", "quantities",
            "missingness", "effects", "design", "lipid_classes",
            "isomer_like_relations"
        )),
        independent_of_package_reader = TRUE,
        support_widening_allowed = FALSE
    )
}

lipidPublicationExecution <- function(dimensions, workload_class) {
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

lipidPublicationPrivacy <- function(route, workload_class) {
    list(
        classification = if (identical(
            workload_class,
            "fixture_correctness"
        ) && identical(route, "lipidsearch")) {
            "public_fixture"
        } else {
            "public_generated"
        },
        private_source_opened = FALSE,
        payload_tracked = FALSE,
        direct_identifiers_present = FALSE,
        cross_omic_values_reused = FALSE
    )
}

lipidPublicationClaimScope <- function(route, profile_id, workload_class) {
    fixture <- identical(workload_class, "fixture_correctness")
    reviewed <- fixture && identical(route, "lipidsearch")
    gc_smoke <- reviewed && profile_id %in% c("gcms", "mixed_lc_gcms")
    performance <- workload_class %in% c(
        "representative",
        "operational_heavy"
    )
    list(
        evidence_class = if (fixture) {
            lipidPublicationFixtureEvidenceClass(route, profile_id)
        } else if (identical(workload_class, "stress")) {
            "stress_characterization"
        } else {
            "public_generated"
        },
        verified_stages = lipidPublicationCapabilities()[[route]]$verified_stages,
        scientific_authority = reviewed,
        performance_authority = performance,
        cross_project_authority = FALSE,
        vendor_authority = reviewed && !gc_smoke,
        gc_vendor_authority = FALSE,
        three_file_workflow_authority = FALSE,
        promotion_authority = FALSE,
        limitations = list(
            "Project-specific evidence without independent real-project sources.",
            "Generated evidence cannot widen the declared route support boundary."
        )
    )
}

lipidPublicationBuildContract <- function(
    route,
    profile_id,
    workload_class,
    fixture_root = "tests/testdata/omics-performance/lipidomics/fixtures"
) {
    paths <- lipidPublicationAuthorityPaths()
    splits <- publicationReadJson(paths$splits)
    assignment <- lipidPublicationSplitAssignmentFor(
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
        lipidPublicationGeneratedDimensions(route, profile_id, workload_class)
    }
    member_count <- dimensions$payload_member_count
    contract <- list(
        schema = "multischolar.omics_publication_lipidomics_workload",
        schema_version = .LIPID_PUBLICATION_VERSION,
        contract_id = paste(
            "lipidomics", route, profile_id, workload_class, "v1", sep = "."
        ),
        owner_ticket_id = .LIPID_PUBLICATION_OWNER,
        status = "frozen",
        workload_id = assignment$assignment_id,
        workload_class = workload_class,
        capability = lipidPublicationExpectedCapability(route),
        assay_profile = list(
            profile_id = profile_id,
            assays = lipidPublicationAssayProfiles()[[profile_id]]$assays,
            payload_mode = if (member_count == 2L) "bundle" else "single",
            member_count = as.integer(member_count)
        ),
        dimensions = dimensions,
        model_profile_id = assignment$generator_parameter_profile_id,
        parameter_authority = lipidPublicationAuthorityBinding(paths$parameters),
        source_authority = lipidPublicationAuthorityBinding(paths$sources),
        split_authority = lipidPublicationAuthorityBinding(paths$splits),
        mapping_authority = lipidPublicationAuthorityBinding(paths$mappings),
        contract_mapping_id = paste(
            "lipidomics.mapping",
            route,
            profile_id,
            sep = "."
        ),
        support_boundary = lipidPublicationAuthorityBinding(paths$boundaries),
        generator = if (fixture) {
            lipidPublicationFixtureGenerator(route, profile_id, fixture_root)
        } else {
            lipidPublicationGeneratedGenerator(route, profile_id)
        },
        rng = lipidPublicationRng(assignment, workload_class),
        truth_contract = lipidPublicationTruthContract(
            route,
            profile_id,
            workload_class
        ),
        execution = lipidPublicationExecution(dimensions, workload_class),
        privacy = lipidPublicationPrivacy(route, workload_class),
        claim_scope = lipidPublicationClaimScope(
            route,
            profile_id,
            workload_class
        ),
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
    lipidPublicationValidateWorkload(contract)
    contract
}

lipidPublicationContractFilename <- function(
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

lipidPublicationFreezeFixtureContracts <- function(output_root) {
    dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
    paths <- character()
    for (route in names(lipidPublicationCapabilities())) {
        for (profile_id in names(lipidPublicationAssayProfiles())) {
            contract <- lipidPublicationBuildContract(
                route,
                profile_id,
                "fixture_correctness"
            )
            path <- file.path(output_root, lipidPublicationContractFilename(
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
