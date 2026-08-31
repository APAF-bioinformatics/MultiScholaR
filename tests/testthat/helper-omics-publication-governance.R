publicationGovernanceRepoPath <- function(...) {
    file.path(
        normalizePath(testthat::test_path("..", ".."), mustWork = TRUE),
        ...
    )
}

publicationGovernanceRecordPath <- function(name) {
    publicationGovernanceRepoPath(
        "tests",
        "testdata",
        "omics-performance",
        name
    )
}

publicationGovernanceRead <- function(name) {
    jsonlite::read_json(
        publicationGovernanceRecordPath(name),
        simplifyVector = FALSE
    )
}

publicationGovernanceReadPath <- function(path) {
    jsonlite::read_json(
        publicationGovernanceRepoPath(path),
        simplifyVector = FALSE
    )
}

publicationGovernanceDigest <- function(path) {
    digest::digest(
        file = publicationGovernanceRepoPath(path),
        algo = "sha256",
        serialize = FALSE
    )
}

publicationGovernanceAbort <- function(message) {
    rlang::abort(
        message,
        class = "multischolar_publication_governance_error"
    )
}

publicationGovernanceScalarString <- function(value) {
    is.character(value) && length(value) == 1L && !is.na(value) && nzchar(value)
}

publicationGovernanceScalarFlag <- function(value) {
    is.logical(value) && length(value) == 1L && !is.na(value)
}

publicationGovernanceScalarNumber <- function(value) {
    is.numeric(value) && length(value) == 1L && is.finite(value)
}

publicationGovernanceDigestValid <- function(value) {
    publicationGovernanceScalarString(value) &&
        grepl("^[0-9a-f]{64}$", value)
}

publicationGovernanceRequireNames <- function(value, expected, label) {
    if (!is.list(value) || !setequal(names(value), expected)) {
        publicationGovernanceAbort(paste(label, "fields differ from contract"))
    }
    invisible(value)
}

publicationGovernanceRequireUnique <- function(values, label) {
    values <- unlist(values, use.names = FALSE)
    if (anyNA(values) || any(!nzchar(values)) || anyDuplicated(values)) {
        publicationGovernanceAbort(paste(label, "must be unique non-empty values"))
    }
    invisible(values)
}

publicationGovernanceHeader <- function(record, schema) {
    if (!is.list(record) || !identical(record$schema, schema) ||
        !identical(record$schema_version, "1.0.0")) {
        publicationGovernanceAbort(paste(schema, "header is invalid"))
    }
    invisible(record)
}

publicationGovernanceCopy <- function(value) {
    jsonlite::fromJSON(
        jsonlite::toJSON(value, auto_unbox = TRUE, null = "null", digits = 17),
        simplifyVector = FALSE
    )
}

publicationGovernanceCapabilities <- function() {
    inventory <- publicationGovernanceReadPath(
        "tests/testdata/omics-capabilities.json"
    )
    capabilities <- unlist(
        lapply(inventory$formats, `[[`, "capabilities"),
        recursive = FALSE
    )
    list(inventory = inventory, capabilities = capabilities)
}

publicationGovernanceCapabilityIds <- function(capabilities) {
    vapply(capabilities, `[[`, character(1), "capability_id")
}

publicationGovernanceValidateReconciliation <- function(record) {
    publicationGovernanceHeader(
        record,
        "multischolar.current_policy_reconciliation"
    )
    publicationGovernanceRequireNames(
        record,
        c(
            "schema", "schema_version", "record_id", "owner_ticket_id",
            "code_revision", "predecessors", "supersession",
            "current_runtime_source", "implemented_decisions",
            "reconciled_current_decision", "rollback", "contradiction_policy"
        ),
        "policy reconciliation"
    )
    predecessor_roles <- vapply(
        record$predecessors,
        `[[`,
        character(1),
        "role"
    )
    if (!identical(
        predecessor_roles,
        c("operational_closeout", "adaptive_auto_policy")
    )) {
        publicationGovernanceAbort("policy predecessor roles are invalid")
    }
    for (predecessor in record$predecessors) {
        if (!publicationGovernanceDigestValid(predecessor$sha256) ||
            !identical(
                publicationGovernanceDigest(predecessor$path),
                predecessor$sha256
            )) {
            publicationGovernanceAbort("policy predecessor digest mismatch")
        }
    }
    for (evidence in record$supersession$later_evidence) {
        if (!identical(
            publicationGovernanceDigest(evidence$path),
            evidence$ticket_sha256
        )) {
            publicationGovernanceAbort("later ticket evidence digest mismatch")
        }
    }
    publicationGovernanceValidateReconciliationRuntime(record)
    publicationGovernanceValidateReconciliationDecisions(record)
    invisible(record)
}

publicationGovernanceValidateReconciliationRuntime <- function(record) {
    source <- record$current_runtime_source
    source_direct <- publicationGovernanceDigestValid(source$sha256) &&
        identical(publicationGovernanceDigest(source$path), source$sha256)
    if (!source_direct &&
        !publicationGovernanceValidateReconciliationTransition(record)) {
        publicationGovernanceAbort("runtime policy source digest mismatch")
    }
    policy <- publicationGovernanceReadPath(
        "tests/testdata/omics-parity/all-omics-auto-policy-v1.json"
    )
    adaptive_ids <- unlist(policy$adaptive_capability_ids, use.names = FALSE)
    decision_ids <- vapply(
        record$implemented_decisions,
        `[[`,
        character(1),
        "capability_id"
    )
    if (!setequal(adaptive_ids, decision_ids) || anyDuplicated(decision_ids)) {
        publicationGovernanceAbort("adaptive capability ownership mismatch")
    }
    if (!identical(
        as.numeric(policy$adaptive_policy$minimum_projected_source_bytes),
        33554432
    )) {
        publicationGovernanceAbort("legacy runtime threshold differs")
    }
    namespace <- asNamespace("MultiScholaR")
    gate_fn <- get(
        record$current_runtime_source$gate_function,
        envir = namespace,
        inherits = FALSE
    )
    if (!identical(
        as.numeric(gate_fn()$minimum_projected_source_bytes),
        33554432
    )) {
        publicationGovernanceAbort("live package gate differs")
    }
    closeout <- publicationGovernanceReadPath(
        "tests/testdata/omics-parity/all-omics-operational-closeout-v1.json"
    )
    if (!identical(closeout$decision, record$predecessors[[1L]]$decision) ||
        !identical(
            closeout$effective_default_backend,
            record$predecessors[[1L]]$effective_default_backend
        ) || length(closeout$promoted_capability_ids)) {
        publicationGovernanceAbort("operational predecessor semantics differ")
    }
    invisible(TRUE)
}

publicationGovernanceValidateReconciliationTransition <- function(record) {
    path <- operationalArtifactRepoPath(
        "inst",
        "extdata",
        "omics-auto-policy-receipts-v2.json"
    )
    if (!file.exists(path)) return(FALSE)
    envelope <- tryCatch(
        workflowPolicyLoadEnvelope(path, strict = TRUE),
        error = function(...) NULL
    )
    if (is.null(envelope)) return(FALSE)
    reconciliation_path <- paste0(
        "tests/testdata/omics-performance/",
        "current-policy-reconciliation-v1.json"
    )
    generated_from_valid <- identical(
        envelope$generated_from$sha256,
        publicationGovernanceDigest(reconciliation_path)
    )
    receipt_ids <- vapply(
        envelope$receipts,
        `[[`,
        character(1),
        "capability_id"
    )
    decision_ids <- vapply(
        record$implemented_decisions,
        `[[`,
        character(1),
        "capability_id"
    )
    receipt_valid <- all(vapply(envelope$receipts, function(receipt) {
        bindings <- Filter(function(binding) {
            identical(binding$binding_id, "legacy_runtime_source")
        }, receipt$authority_bindings)
        length(bindings) == 1L &&
            identical(bindings[[1L]]$path_or_uri, record$current_runtime_source$path) &&
            identical(bindings[[1L]]$sha256, record$current_runtime_source$sha256) &&
            identical(receipt$receipt_kind, "reconciled_current") &&
            identical(receipt$decision, "engineering_current_baseline") &&
            !isTRUE(receipt$claim_scope$publication_authority)
    }, logical(1)))
    isTRUE(generated_from_valid) && setequal(receipt_ids, decision_ids) &&
        !anyDuplicated(receipt_ids) && receipt_valid
}

publicationGovernanceValidateReconciliationDecisions <- function(record) {
    namespace <- asNamespace("MultiScholaR")
    descriptor_fns <- c(
        "metabolomics.custom.metabolite.standard.v1" =
            "artifactMetabolomicsWorkflowDescriptor",
        "lipidomics.lipidsearch.lipid.standard.v1" =
            "artifactLipidomicsWorkflowDescriptor"
    )
    closeout <- publicationGovernanceReadPath(
        "tests/testdata/omics-parity/all-omics-operational-closeout-v1.json"
    )
    predecessor_digests <- stats::setNames(
        vapply(
            closeout$descriptor_canaries,
            `[[`,
            character(1),
            "descriptor_digest"
        ),
        vapply(
            closeout$descriptor_canaries,
            `[[`,
            character(1),
            "capability_id"
        )
    )
    for (decision in record$implemented_decisions) {
        dynamic <- decision$dynamic_runtime_decision
        if (!identical(decision$static_certification$status, "evict") ||
            isTRUE(decision$static_certification$auto_eligible) ||
            !identical(dynamic$platform, "Linux") ||
            !identical(as.numeric(dynamic$threshold_bytes), 33554432) ||
            !isTRUE(dynamic$legacy_size_measure_is_heuristic) ||
            isTRUE(decision$engineering_evidence$publication_authority)) {
            publicationGovernanceAbort("implemented decision is contradictory")
        }
        descriptor_fn <- get(
            descriptor_fns[[decision$capability_id]],
            envir = namespace,
            inherits = FALSE
        )
        if (!identical(descriptor_fn()$descriptor_digest, decision$descriptor_digest)) {
            publicationGovernanceAbort("live descriptor digest mismatch")
        }
        if (!identical(
            predecessor_digests[[decision$capability_id]],
            decision$predecessor_descriptor_digest
        )) {
            publicationGovernanceAbort("predecessor descriptor digest mismatch")
        }
        for (binding in list(decision$descriptor_source, decision$scale_workload)) {
            if (!identical(publicationGovernanceDigest(binding$path), binding$sha256)) {
                publicationGovernanceAbort("decision binding digest mismatch")
            }
        }
    }
    if (isTRUE(record$reconciled_current_decision$publication_decision_supported)) {
        publicationGovernanceAbort("engineering reconciliation cannot publish")
    }
    invisible(TRUE)
}

publicationGovernanceValidateProtocol <- function(record) {
    publicationGovernanceHeader(record, "multischolar.omics_publication_protocol")
    publicationGovernanceRequireNames(
        record,
        c(
            "schema", "schema_version", "protocol_id", "owner_ticket_id",
            "status", "authorities", "comparators", "governance_records",
            "claim_classes", "workload_classes", "independence", "replication",
            "cache_policy", "host_policy", "primary_endpoints",
            "secondary_endpoints", "statistics", "promotion_gates", "privacy",
            "terminal_outcomes", "amendment_policy"
        ),
        "publication protocol"
    )
    revisions <- vapply(record$comparators[1:2], `[[`, character(1), "revision")
    if (!identical(
        revisions,
        c(
            "c7c12851ea4d9e96f91df0e1f9e7b91c2eb51017",
            "56725e90c5fda97835775f5f8f57c02703b53120"
        )
    )) {
        publicationGovernanceAbort("comparator revisions differ")
    }
    publicationGovernanceValidateProtocolDependencies(record)
    publicationGovernanceValidateProtocolDesign(record)
    invisible(record)
}

publicationGovernanceValidateProtocolDependencies <- function(record) {
    for (authority in record$authorities) {
        if (!file.exists(publicationGovernanceRepoPath(authority$path))) {
            publicationGovernanceAbort("protocol authority is missing")
        }
        if (!is.null(authority$sha256) &&
            !identical(publicationGovernanceDigest(authority$path), authority$sha256)) {
            publicationGovernanceAbort("protocol authority digest mismatch")
        }
    }
    paths <- vapply(record$governance_records, `[[`, character(1), "path")
    publicationGovernanceRequireUnique(paths, "governance record paths")
    if (!all(file.exists(vapply(paths, publicationGovernanceRepoPath, character(1))))) {
        publicationGovernanceAbort("protocol governance record is missing")
    }
    invisible(TRUE)
}

publicationGovernanceValidateManifest <- function(record) {
    publicationGovernanceHeader(
        record,
        "multischolar.omics_publication_governance_manifest"
    )
    publicationGovernanceRequireNames(
        record,
        c(
            "schema", "schema_version", "manifest_id", "owner_ticket_id",
            "records", "immutability"
        ),
        "publication governance manifest"
    )
    paths <- vapply(record$records, `[[`, character(1), "path")
    publicationGovernanceRequireUnique(paths, "governance manifest paths")
    for (entry in record$records) {
        if (!publicationGovernanceDigestValid(entry$sha256) ||
            !identical(publicationGovernanceDigest(entry$path), entry$sha256)) {
            publicationGovernanceAbort("governance manifest digest mismatch")
        }
    }
    if (!isTRUE(record$immutability$additive_successor_required) ||
        isTRUE(record$immutability$mutate_version_in_place) ||
        !isTRUE(record$immutability$predecessor_replay_required)) {
        publicationGovernanceAbort("governance manifest is mutable")
    }
    invisible(record)
}

publicationGovernanceValidateProtocolDesign <- function(record) {
    independence <- record$independence
    replication <- record$replication
    if (independence$minimum_projects_per_cross_project_claim < 3 ||
        independence$minimum_total_primary_pairs < 30 ||
        independence$minimum_pairs_per_project < 10 ||
        isTRUE(independence$generated_workload_counts_as_real_project) ||
        replication$warmup_runs_per_backend_project_estimand != 2 ||
        isTRUE(replication$optional_stopping) ||
        isTRUE(record$amendment_policy$candidate_result_may_change_protocol)) {
        publicationGovernanceAbort("publication experimental design is weakened")
    }
    if (!is.null(replication$primary_pair_count) ||
        replication$minimum_primary_pair_count != 30 ||
        replication$maximum_primary_pair_count != 60 ||
        !identical(
            replication$precision_selection_owner_ticket_id,
            "OMICS-ART-063"
        )) {
        publicationGovernanceAbort("pair-count precision ownership differs")
    }
    required_classes <- c(
        "fixture_correctness", "representative", "operational_heavy", "stress"
    )
    classes <- vapply(record$workload_classes, `[[`, character(1), "workload_class")
    if (!identical(classes, required_classes)) {
        publicationGovernanceAbort("workload classes differ")
    }
    if (!all(vapply(record$promotion_gates, publicationGovernanceScalarNumber, logical(1)))) {
        non_numeric <- c("all_scientific_and_lifecycle_gates_required")
        numeric_gates <- record$promotion_gates[setdiff(
            names(record$promotion_gates),
            non_numeric
        )]
        if (!all(vapply(numeric_gates, publicationGovernanceScalarNumber, logical(1))) ||
            !isTRUE(record$promotion_gates$all_scientific_and_lifecycle_gates_required)) {
            publicationGovernanceAbort("promotion gates are invalid")
        }
    }
    invisible(TRUE)
}

publicationGovernanceValidateCoverage <- function(record) {
    publicationGovernanceHeader(record, "multischolar.omics_publication_coverage")
    publicationGovernanceRequireNames(
        record,
        c(
            "schema", "schema_version", "coverage_id", "owner_ticket_id",
            "capability_inventory", "full_workflows", "reader_boundaries",
            "excluded_capabilities", "excluded_formats", "auxiliary_surfaces",
            "placeholder_surfaces", "counts", "unknown_or_unowned_policy"
        ),
        "publication coverage"
    )
    source <- publicationGovernanceCapabilities()
    inventory <- source$inventory
    capability_ids <- publicationGovernanceCapabilityIds(source$capabilities)
    covered_ids <- c(
        vapply(record$full_workflows, `[[`, character(1), "capability_id"),
        vapply(record$reader_boundaries, `[[`, character(1), "capability_id"),
        vapply(record$excluded_capabilities, `[[`, character(1), "capability_id")
    )
    if (!setequal(capability_ids, covered_ids) || anyDuplicated(covered_ids)) {
        publicationGovernanceAbort("capability coverage is incomplete or duplicate")
    }
    publicationGovernanceValidateCoverageCapabilities(record, source$capabilities)
    publicationGovernanceValidateCoverageFormats(record, inventory)
    publicationGovernanceValidateCoverageOwners(record)
    invisible(record)
}

publicationGovernanceValidateCoverageCapabilities <- function(
    record,
    capabilities
) {
    by_id <- stats::setNames(
        capabilities,
        publicationGovernanceCapabilityIds(capabilities)
    )
    for (entry in record$full_workflows) {
        capability <- by_id[[entry$capability_id]]
        if (is.null(capability) || !isTRUE(capability$complete_workflow) ||
            !identical(entry$support_status, "scientifically_supported") ||
            !setequal(
                unlist(entry$e2e_lanes, use.names = FALSE),
                unlist(capability$e2e_lanes, use.names = FALSE)
            )) {
            publicationGovernanceAbort("full workflow coverage semantics differ")
        }
    }
    for (entry in record$reader_boundaries) {
        capability <- by_id[[entry$capability_id]]
        if (is.null(capability) || isTRUE(capability$complete_workflow) ||
            !identical(entry$support_status, "reader_characterized") ||
            isTRUE(entry$full_workflow_promotion_authority)) {
            publicationGovernanceAbort("reader boundary coverage semantics differ")
        }
    }
    for (entry in record$excluded_capabilities) {
        capability <- by_id[[entry$capability_id]]
        if (is.null(capability) || isTRUE(capability$complete_workflow) ||
            !identical(entry$support_status, "advertised_unverified") ||
            isTRUE(entry$positive_benchmark_authority)) {
            publicationGovernanceAbort("excluded capability semantics differ")
        }
    }
    invisible(TRUE)
}

publicationGovernanceValidateCoverageFormats <- function(record, inventory) {
    detection_ids <- vapply(
        Filter(
            \(format) identical(format$support_status, "detection_only"),
            inventory$formats
        ),
        `[[`,
        character(1),
        "format_id"
    )
    excluded_ids <- vapply(record$excluded_formats, `[[`, character(1), "format_id")
    if (!setequal(detection_ids, excluded_ids) || anyDuplicated(excluded_ids)) {
        publicationGovernanceAbort("detection-only format coverage differs")
    }
    surface_ids <- vapply(
        inventory$non_workflow_surfaces,
        `[[`,
        character(1),
        "surface_id"
    )
    covered_surfaces <- c(
        vapply(record$auxiliary_surfaces, `[[`, character(1), "surface_id"),
        vapply(record$placeholder_surfaces, `[[`, character(1), "surface_id")
    )
    if (!setequal(surface_ids, covered_surfaces) || anyDuplicated(covered_surfaces)) {
        publicationGovernanceAbort("non-workflow surface coverage differs")
    }
    invisible(TRUE)
}

publicationGovernanceValidateCoverageOwners <- function(record) {
    tickets <- jsonlite::read_json(
        publicationGovernanceRepoPath("tickets/openai/tickets-index.json"),
        simplifyVector = FALSE
    )$tickets
    ticket_ids <- vapply(tickets, `[[`, character(1), "id")
    owners <- unlist(lapply(record$full_workflows, \(entry) {
        c(
            entry$workload_owner_ticket_id,
            unlist(entry$runtime_owner_ticket_ids, use.names = FALSE),
            entry$decision_owner_ticket_id
        )
    }), use.names = FALSE)
    if (!all(owners %in% ticket_ids) ||
        any(vapply(record$full_workflows, \(entry) {
            entry$minimum_independent_projects < 3
        }, logical(1)))) {
        publicationGovernanceAbort("coverage owner or project rule is invalid")
    }
    invisible(TRUE)
}

publicationGovernanceValidateProjects <- function(record, coverage) {
    publicationGovernanceHeader(record, "multischolar.omics_publication_projects")
    publicationGovernanceRequireNames(
        record,
        c(
            "schema", "schema_version", "projects_id", "owner_ticket_id",
            "minimum_real_project_authorities_per_cross_project_claim",
            "generated_payload_counts_as_real_project",
            "aggregate_calibrated_generated_payload_counts_as_real_project",
            "source_receipt_required_fields", "full_workflow_claims",
            "reader_boundary_claims", "auxiliary_claims", "fixture_authority",
            "successor_policy"
        ),
        "publication projects"
    )
    full_ids <- vapply(coverage$full_workflows, `[[`, character(1), "capability_id")
    project_ids <- vapply(
        record$full_workflow_claims,
        `[[`,
        character(1),
        "capability_id"
    )
    if (!identical(full_ids, project_ids) ||
        isTRUE(record$generated_payload_counts_as_real_project)) {
        publicationGovernanceAbort("project capability or synthetic rule differs")
    }
    for (claim in record$full_workflow_claims) {
        if (claim$verified_real_project_count < claim$required_real_project_count &&
            (isTRUE(claim$promotion_eligible) ||
                !grepl("insufficient", claim$cross_project_claim_status))) {
            publicationGovernanceAbort("insufficient projects cannot promote")
        }
        kinds <- vapply(
            claim$current_real_project_authorities,
            `[[`,
            character(1),
            "source_kind"
        )
        if (any(grepl("generated", kinds))) {
            publicationGovernanceAbort("generated source occupies real project slot")
        }
    }
    if (isTRUE(record$successor_policy$candidate_access_before_successor)) {
        publicationGovernanceAbort("candidate cannot precede project successor")
    }
    invisible(record)
}

publicationGovernanceValidateRoles <- function(record) {
    publicationGovernanceHeader(record, "multischolar.omics_publication_roles")
    publicationGovernanceRequireNames(
        record,
        c(
            "schema", "schema_version", "roles_id", "owner_ticket_id",
            "roles", "separation_rules", "handoffs", "blinding", "readiness"
        ),
        "publication roles"
    )
    role_ids <- vapply(record$roles, `[[`, character(1), "role_id")
    expected <- c(
        "protocol_owner", "source_data_owner", "implementation_owner",
        "benchmark_operator", "analysis_owner", "independent_reviewer"
    )
    if (!identical(role_ids, expected) || anyDuplicated(role_ids) ||
        isTRUE(record$separation_rules$self_approval_allowed)) {
        publicationGovernanceAbort("role separation contract differs")
    }
    assigned <- vapply(
        record$roles,
        \(role) publicationGovernanceScalarString(role$principal_id),
        logical(1)
    )
    if (!all(assigned) &&
        (isTRUE(record$readiness$runtime_implementation_authorized) ||
            isTRUE(record$readiness$campaign_execution_authorized))) {
        publicationGovernanceAbort("unassigned roles cannot authorize work")
    }
    if (all(assigned)) {
        principals <- stats::setNames(
            vapply(record$roles, `[[`, character(1), "principal_id"),
            role_ids
        )
        if (principals[["independent_reviewer"]] %in%
            principals[c("implementation_owner", "analysis_owner")]) {
            publicationGovernanceAbort("independent reviewer conflicts")
        }
    }
    invisible(record)
}

publicationGovernanceValidateEstimands <- function(record, coverage) {
    publicationGovernanceHeader(record, "multischolar.omics_publication_estimands")
    publicationGovernanceRequireNames(
        record,
        c(
            "schema", "schema_version", "estimands_id", "owner_ticket_id",
            "primary_import_work_unit", "effect_classes", "phase_definitions",
            "capability_phase_assignments", "concrete_id_rule", "consumers"
        ),
        "publication estimands"
    )
    phase_ids <- vapply(record$phase_definitions, `[[`, character(1), "phase_id")
    expected_phases <- c(
        "import_and_settle", "scientific_stage", "complete_workflow",
        "bounded_query", "closed_project_resume", "auxiliary_api_call"
    )
    if (!identical(phase_ids, expected_phases) || anyDuplicated(phase_ids)) {
        publicationGovernanceAbort("estimand phase set differs")
    }
    for (phase in record$phase_definitions) {
        required <- c(
            "phase_id", "start_marker", "stop_marker", "included_operations",
            "excluded_operations", "artifact_precondition", "cache_precondition",
            "primary_work_unit_id", "secondary_work_unit_ids",
            "allowed_effect_classes", "terminal_failures"
        )
        publicationGovernanceRequireNames(phase, required, phase$phase_id)
        if (!publicationGovernanceScalarString(phase$start_marker) ||
            !publicationGovernanceScalarString(phase$stop_marker)) {
            publicationGovernanceAbort("estimand markers are invalid")
        }
    }
    publicationGovernanceValidateEstimandAssignments(record, coverage, phase_ids)
    invisible(record)
}

publicationGovernanceValidateEstimandAssignments <- function(
    record,
    coverage,
    phase_ids
) {
    expected_subjects <- c(
        vapply(coverage$full_workflows, `[[`, character(1), "capability_id"),
        vapply(coverage$reader_boundaries, `[[`, character(1), "capability_id"),
        vapply(coverage$auxiliary_surfaces, `[[`, character(1), "surface_id")
    )
    observed <- vapply(
        record$capability_phase_assignments,
        `[[`,
        character(1),
        "subject_id"
    )
    if (!setequal(expected_subjects, observed) || anyDuplicated(observed)) {
        publicationGovernanceAbort("estimand subject coverage differs")
    }
    for (assignment in record$capability_phase_assignments) {
        assigned_phases <- unlist(assignment$phase_ids, use.names = FALSE)
        if (!all(assigned_phases %in% phase_ids) ||
            (isTRUE(assignment$complete_workflow) &&
                !"complete_workflow" %in% assigned_phases)) {
            publicationGovernanceAbort("estimand assignment is invalid")
        }
    }
    invisible(TRUE)
}

publicationGovernanceValidateSplits <- function(record, projects) {
    publicationGovernanceHeader(record, "multischolar.omics_publication_splits")
    publicationGovernanceRequireNames(
        record,
        c(
            "schema", "schema_version", "splits_id", "owner_ticket_id",
            "split_units", "disjointness_rules", "generated_seed_families",
            "assignments", "readiness"
        ),
        "publication splits"
    )
    families <- record$generated_seed_families
    ranges <- lapply(families, \(family) {
        seq.int(family$minimum_seed, family$maximum_seed)
    })
    if (length(intersect(ranges$pilot, ranges$holdout)) ||
        length(intersect(ranges$pilot, ranges$stress)) ||
        length(intersect(ranges$holdout, ranges$stress))) {
        publicationGovernanceAbort("generated seed families overlap")
    }
    for (assignment in record$assignments) {
        pilot_projects <- unlist(assignment$pilot_project_ids, use.names = FALSE)
        holdout_projects <- unlist(assignment$holdout_project_ids, use.names = FALSE)
        pilot_sources <- unlist(assignment$pilot_source_ids, use.names = FALSE)
        holdout_sources <- unlist(assignment$holdout_source_ids, use.names = FALSE)
        if (length(intersect(pilot_projects, holdout_projects)) ||
            length(intersect(pilot_sources, holdout_sources))) {
            publicationGovernanceAbort("pilot and holdout assignments overlap")
        }
        if ((!length(pilot_projects) || !length(holdout_projects)) &&
            isTRUE(assignment$promotion_eligible)) {
            publicationGovernanceAbort("empty split cannot promote")
        }
    }
    if (!isTRUE(record$readiness$successor_required_before_candidate) ||
        isTRUE(record$readiness$candidate_access_allowed)) {
        publicationGovernanceAbort("split readiness is unsafe")
    }
    invisible(record)
}

publicationGovernanceValidateThresholdGrid <- function(record) {
    publicationGovernanceHeader(
        record,
        "multischolar.omics_publication_threshold_grid"
    )
    publicationGovernanceRequireNames(
        record,
        c(
            "schema", "schema_version", "threshold_grid_id", "owner_ticket_id",
            "size_measure", "threshold_bytes", "selection_rule",
            "heldout_validation", "terminal_outcomes", "promotion_policy",
            "legacy_current_policy"
        ),
        "publication threshold grid"
    )
    measure <- record$size_measure
    if (!identical(measure$measure_id, "total_uncompressed_input_bytes_v1") ||
        !isTRUE(measure$backend_invariant) ||
        isTRUE(measure$uses_materialized_r_payload) ||
        isTRUE(measure$uses_object_size_multiplier) ||
        isTRUE(measure$uses_fitted_expansion_factor)) {
        publicationGovernanceAbort("threshold size measure is heuristic")
    }
    grid <- as.numeric(unlist(record$threshold_bytes, use.names = FALSE))
    if (any(!is.finite(grid)) || any(grid <= 0) ||
        is.unsorted(grid, strictly = TRUE) || anyDuplicated(grid)) {
        publicationGovernanceAbort("threshold grid is invalid")
    }
    if (!isTRUE(record$selection_rule$monotonic_persistence_required) ||
        isTRUE(record$heldout_validation$confirmatory_results_may_reselect_threshold) ||
        !isTRUE(record$promotion_policy$stable_threshold_validated_required)) {
        publicationGovernanceAbort("threshold decision policy is unsafe")
    }
    invisible(record)
}

publicationGovernanceValidateRetry <- function(record) {
    publicationGovernanceHeader(
        record,
        "multischolar.omics_publication_retry_policy"
    )
    publicationGovernanceRequireNames(
        record,
        c(
            "schema", "schema_version", "retry_policy_id", "owner_ticket_id",
            "attempt_contract", "retryable_failures",
            "non_retryable_candidate_failures",
            "non_retryable_campaign_failures", "slot_outcomes", "campaign_rule"
        ),
        "publication retry policy"
    )
    attempts <- record$attempt_contract
    if (attempts$initial_attempts_per_schedule_slot != 1 ||
        attempts$maximum_retries_per_schedule_slot != 1 ||
        attempts$maximum_attempts_per_schedule_slot != 2 ||
        isTRUE(attempts$overwrite_attempt_record) ||
        !isTRUE(attempts$append_only_ledger)) {
        publicationGovernanceAbort("retry attempt bounds differ")
    }
    if (any(vapply(record$retryable_failures, \(failure) {
        isTRUE(failure$worker_pid_existed) || isTRUE(failure$candidate_work_executed)
    }, logical(1)))) {
        publicationGovernanceAbort("candidate work cannot be retried")
    }
    if (!isTRUE(record$campaign_rule$candidate_failure_counts_against_decision) ||
        isTRUE(record$campaign_rule$failed_slot_replaced_by_new_slot) ||
        isTRUE(record$campaign_rule$post_hoc_attempt_extension)) {
        publicationGovernanceAbort("retry campaign rule is selective")
    }
    invisible(record)
}

publicationGovernanceValidateCampaignMatrix <- function(record) {
    publicationGovernanceHeader(
        record,
        "multischolar.omics_publication_campaign_matrix"
    )
    publicationGovernanceRequireNames(
        record,
        c(
            "schema", "schema_version", "campaign_matrix_id", "owner_ticket_id",
            "status", "pair_design", "subject_groups", "secondary_host_reserve",
            "setup_and_parity_reserve", "totals", "current_readiness"
        ),
        "publication campaign matrix"
    )
    group_max <- sum(vapply(
        record$subject_groups,
        \(group) as.numeric(group$maximum_group_process_launches),
        numeric(1)
    ))
    totals <- record$totals
    computed <- group_max +
        as.numeric(record$secondary_host_reserve$maximum_process_launches) +
        as.numeric(record$setup_and_parity_reserve$maximum_process_launches)
    if (!identical(computed, as.numeric(totals$computed_maximum_process_launches)) ||
        totals$governed_maximum_process_launches < computed ||
        isTRUE(record$current_readiness$campaign_execution_authorized)) {
        publicationGovernanceAbort("campaign matrix totals or readiness differ")
    }
    for (group in record$subject_groups) {
        pair_slots <- sum(vapply(
            group$paired_strata_per_subject,
            \(stratum) as.numeric(stratum$maximum_pair_slots),
            numeric(1)
        ))
        if (!identical(
            pair_slots,
            as.numeric(group$maximum_pair_slots_per_subject)
        ) || !identical(
            as.numeric(group$maximum_measured_process_slots_per_subject),
            2 * pair_slots
        ) ||
            !identical(
                group$maximum_group_process_launches,
                group$subject_count * group$maximum_process_launches_per_subject
            )) {
            publicationGovernanceAbort("campaign group arithmetic differs")
        }
    }
    invisible(record)
}

publicationGovernanceValidateBudget <- function(record, matrix, retry) {
    publicationGovernanceHeader(
        record,
        "multischolar.omics_publication_campaign_budget"
    )
    publicationGovernanceRequireNames(
        record,
        c(
            "schema", "schema_version", "campaign_budget_id", "owner_ticket_id",
            "campaign_matrix", "hard_limits", "host_preflight_required",
            "mid_run_abort_conditions", "resume_contract",
            "retention_and_cleanup", "current_status"
        ),
        "publication campaign budget"
    )
    hard <- record$hard_limits
    if (hard$maximum_process_launches <
            matrix$totals$computed_maximum_process_launches ||
        hard$maximum_process_launches >
            matrix$totals$governed_maximum_process_launches ||
        hard$maximum_attempts_per_schedule_slot !=
            retry$attempt_contract$maximum_attempts_per_schedule_slot ||
        any(!vapply(hard, publicationGovernanceScalarNumber, logical(1))) ||
        isTRUE(record$current_status$execution_authorized)) {
        publicationGovernanceAbort("campaign budget is invalid or executable")
    }
    if (!all(c(
        "owned_cgroup_v2_available",
        "free_memory_safety_margin_satisfied",
        "free_disk_safety_margin_satisfied",
        "no_thermal_throttling"
    ) %in% unlist(record$host_preflight_required, use.names = FALSE))) {
        publicationGovernanceAbort("campaign host preflight is incomplete")
    }
    invisible(record)
}

publicationGovernanceValidateReceiptSchema <- function(record) {
    required_top <- c(
        "$schema", "$id", "title", "description", "type",
        "additionalProperties", "required", "properties", "$defs", "allOf"
    )
    publicationGovernanceRequireNames(record, required_top, "policy receipt schema")
    required <- unlist(record$required, use.names = FALSE)
    properties <- names(record$properties)
    if (!all(required %in% properties) ||
        isTRUE(record$additionalProperties) ||
        !identical(
            record$properties$size_measure$properties$measure_id$const,
            "total_uncompressed_input_bytes_v1"
        ) ||
        !isTRUE(record$properties$size_measure$properties$exact$const) ||
        !isTRUE(
            record$properties$size_measure$properties$available_before_full_parse$const
        )) {
        publicationGovernanceAbort("policy receipt schema is weakened")
    }
    decision_values <- unlist(
        record$properties$decision$enum,
        use.names = FALSE
    )
    required_decisions <- c(
        "reconciled_current", "proposed_pilot",
        "confirmatory_decision", "final_installed"
    )
    kinds <- unlist(record$properties$receipt_kind$enum, use.names = FALSE)
    if (!setequal(kinds, required_decisions) ||
        !all(c("confirmed_promote", "confirmed_no_promotion") %in%
            decision_values)) {
        publicationGovernanceAbort("policy receipt decision states differ")
    }
    invisible(record)
}

publicationGovernanceValidateReceipt <- function(receipt, schema) {
    required <- unlist(schema$required, use.names = FALSE)
    if (!is.list(receipt) || !all(required %in% names(receipt)) ||
        any(!names(receipt) %in% names(schema$properties))) {
        publicationGovernanceAbort("policy receipt fields differ")
    }
    if (!identical(receipt$schema, "multischolar.omics_auto_policy_receipt") ||
        !identical(receipt$schema_version, "1.0.0") ||
        !receipt$receipt_kind %in% unlist(
            schema$properties$receipt_kind$enum,
            use.names = FALSE
        ) ||
        !receipt$decision %in% unlist(
            schema$properties$decision$enum,
            use.names = FALSE
        )) {
        publicationGovernanceAbort("policy receipt identity differs")
    }
    measure <- receipt$size_measure
    if (!identical(measure$measure_id, "total_uncompressed_input_bytes_v1") ||
        !identical(measure$unit, "byte") || !isTRUE(measure$exact) ||
        !isTRUE(measure$available_before_full_parse)) {
        publicationGovernanceAbort("policy receipt size measure is invalid")
    }
    publicationGovernanceValidateReceiptBindings(receipt)
    publicationGovernanceValidateReceiptDecision(receipt)
    invisible(receipt)
}

publicationGovernanceValidateReceiptBindings <- function(receipt) {
    bindings <- c(receipt$authority_bindings, receipt$evidence_bindings)
    for (binding in bindings) {
        if (!publicationGovernanceScalarString(binding$binding_id) ||
            !publicationGovernanceScalarString(binding$path_or_uri) ||
            !publicationGovernanceDigestValid(binding$sha256)) {
            publicationGovernanceAbort("policy receipt binding is invalid")
        }
    }
    if (!length(receipt$role_bindings) || any(vapply(
        receipt$role_bindings,
        \(binding) {
            !publicationGovernanceScalarString(binding$principal_id) ||
                !publicationGovernanceDigestValid(binding$handoff_digest)
        },
        logical(1)
    ))) {
        publicationGovernanceAbort("policy receipt role binding is invalid")
    }
    if (!publicationGovernanceDigestValid(receipt$receipt_digest)) {
        publicationGovernanceAbort("policy receipt digest is invalid")
    }
    invisible(TRUE)
}

publicationGovernanceValidateReceiptDecision <- function(receipt) {
    promoted <- receipt$decision %in% c("confirmed_promote", "final_installed")
    if (promoted) {
        required <- c(
            "project_ids", "pilot_split_id", "holdout_split_id",
            "threshold_grid_id", "estimand_ids", "candidate_revision",
            "gate_results"
        )
        if (!all(required %in% names(receipt)) ||
            !publicationGovernanceScalarNumber(receipt$threshold_bytes) ||
            !identical(receipt$at_or_above_threshold$backend, "artifact") ||
            !identical(receipt$at_or_above_threshold$rollout, "evict") ||
            !isTRUE(receipt$claim_scope$publication_authority) ||
            length(receipt$project_ids) < 3L) {
            publicationGovernanceAbort("promoted policy receipt is incomplete")
        }
    }
    if (identical(receipt$decision, "confirmed_no_promotion") &&
        (!publicationGovernanceScalarString(receipt$terminal_reason) ||
            !identical(receipt$at_or_above_threshold$backend, "memory") ||
            !identical(receipt$at_or_above_threshold$rollout, "none"))) {
        publicationGovernanceAbort("no-promotion receipt is invalid")
    }
    invisible(TRUE)
}
