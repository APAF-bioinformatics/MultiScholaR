metabPublicationHeavyDimensions <- function(
    route,
    profile_id,
    aggregate_feature_count,
    sample_count
) {
    counts <- metabPublicationAssayCounts(profile_id, aggregate_feature_count)
    member_count <- if (identical(route, "msdial") &&
        identical(profile_id, "mixed")) 3L else 1L
    list(
        aggregate_feature_count = as.integer(aggregate_feature_count),
        assay_feature_counts = counts,
        sample_count = as.integer(sample_count),
        assay_count = as.integer(length(counts)),
        quantity_count = as.numeric(aggregate_feature_count * sample_count),
        payload_member_count = member_count
    )
}

metabPublicationBuildHeavyContract <- function(
    route,
    profile_id,
    aggregate_feature_count,
    sample_count,
    qualification_id
) {
    if (!publicationScalarString(qualification_id)) {
        metabPublicationAbort("heavy qualification identity is invalid")
    }
    contract <- metabPublicationBuildContract(
        route,
        profile_id,
        "operational_heavy"
    )
    dimensions <- metabPublicationHeavyDimensions(
        route,
        profile_id,
        aggregate_feature_count,
        sample_count
    )
    contract$dimensions <- dimensions
    contract$assay_profile$payload_mode <- if (
        dimensions$payload_member_count == 3L
    ) "bundle" else "single"
    contract$assay_profile$member_count <- dimensions$payload_member_count
    contract$model_profile_id <- paste(
        contract$model_profile_id,
        qualification_id,
        sep = "."
    )
    contract$execution <- metabPublicationExecution(
        dimensions,
        "operational_heavy"
    )
    manifest_path <- "tools/profiling/omics_publication_metabolomics_manifest.R"
    contract$generator$source_bindings[[
        length(contract$generator$source_bindings) + 1L
    ]] <- metabPublicationAuthorityBinding(manifest_path)
    contract$expected_digests <- list(
        payload_set_sha256 = strrep("0", 64L),
        truth_sha256 = strrep("0", 64L)
    )
    metabPublicationValidateWorkload(contract)
    contract
}

metabPublicationManifestBindingNames <- function() {
    c(
        "sources", "source_metadata_search", "splits", "parameters",
        "support_boundaries",
        "mapping_authorities", "msdial_bundles", "exclusions", "negative",
        "governance_successor", "fixture_import_parity",
        "representative_determinism", "representative_import_parity",
        "stress_determinism", "stress_import_parity"
    )
}

metabPublicationValidateSourceMetadataSearch <- function(record) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "owner_ticket_id", "status", "searched_at",
        "repositories", "queries", "candidate_count", "admitted_source_count",
        "candidates", "source_bytes_downloaded", "schema_bytes_inspected",
        "license_use_cleared", "independence_promoted", "decision",
        "limitations", "promotion_authority", "publication_authority"
    ), "Metabolomics source metadata search")
    ids <- vapply(record$candidates, `[[`, character(1), "project_id")
    valid <- identical(
        record$schema,
        "multischolar.omics_publication_metabolomics_metadata_search"
    ) && identical(record$schema_version, "1.0.0") &&
        identical(record$owner_ticket_id, .METAB_PUBLICATION_OWNER) &&
        identical(record$status, "metadata_only_nonadmission") &&
        identical(record$candidate_count, as.integer(length(record$candidates))) &&
        record$candidate_count >= 3L && !anyDuplicated(ids) &&
        identical(record$admitted_source_count, 0L) &&
        !isTRUE(record$source_bytes_downloaded) &&
        !isTRUE(record$schema_bytes_inspected) &&
        !isTRUE(record$license_use_cleared) &&
        !isTRUE(record$independence_promoted) &&
        identical(record$decision, "project_specific_nonpromotional") &&
        length(record$limitations) > 0L &&
        !isTRUE(record$promotion_authority) &&
        !isTRUE(record$publication_authority)
    if (!valid) metabPublicationAbort("source metadata search differs")
    invisible(record)
}

metabPublicationValidateManifestWorkload <- function(entry) {
    publicationRequireNames(entry, c(
        "workload_id", "route", "profile_id", "workload_class", "contract",
        "dimensions", "payload_set_sha256", "truth_sha256"
    ), "Metabolomics manifest workload")
    metabPublicationValidateBinding(entry$contract, "Manifest contract")
    contract <- metabPublicationReadContract(entry$contract$path)
    valid <- identical(entry$workload_id, contract$workload_id) &&
        identical(entry$route, contract$capability$input_format) &&
        identical(entry$profile_id, contract$assay_profile$profile_id) &&
        identical(entry$workload_class, contract$workload_class) &&
        identical(entry$dimensions, contract$dimensions) &&
        identical(
            entry$payload_set_sha256,
            contract$expected_digests$payload_set_sha256
        ) && identical(
            entry$truth_sha256,
            contract$expected_digests$truth_sha256
        )
    metabPublicationRequireDigest(entry$payload_set_sha256, "Manifest payload")
    metabPublicationRequireDigest(entry$truth_sha256, "Manifest truth")
    if (!valid) metabPublicationAbort("manifest workload differs")
    contract
}

metabPublicationValidateQualifications <- function(values) {
    routes <- vapply(values, `[[`, character(1), "route")
    valid <- length(values) == 2L && setequal(
        routes,
        names(metabPublicationCapabilities())
    ) && !anyDuplicated(routes) && all(vapply(values, function(value) {
        setequal(names(value), c(
            "route", "qualification_id", "pilot_receipt",
            "aggregate_feature_count", "sample_count", "qualified",
            "candidate_loaded", "performance_authority"
        )) && publicationScalarString(value$qualification_id) &&
            isTRUE(value$qualified) && !isTRUE(value$candidate_loaded) &&
            !isTRUE(value$performance_authority) &&
            value$aggregate_feature_count >= 80000L && value$sample_count >= 48L &&
            tryCatch({
                metabPublicationValidateBinding(
                    value$pilot_receipt,
                    "Pilot qualification"
                )
                TRUE
            }, error = function(error) FALSE)
    }, logical(1)))
    if (!valid) metabPublicationAbort("manifest qualification differs")
    invisible(values)
}

metabPublicationValidateFinalManifest <- function(record) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "manifest_id", "owner_ticket_id", "status",
        "bindings", "pilot_qualifications", "workloads", "source_claim_scope",
        "candidate_access_allowed", "promotion_authority",
        "publication_authority"
    ), "Metabolomics final manifest")
    publicationRequireNames(
        record$bindings,
        metabPublicationManifestBindingNames(),
        "Metabolomics manifest bindings"
    )
    lapply(record$bindings, function(binding) {
        metabPublicationValidateBinding(binding, "Manifest evidence")
    })
    metabPublicationValidateSourceMetadataSearch(publicationReadJson(
        record$bindings$source_metadata_search$path
    ))
    metabPublicationValidateQualifications(record$pilot_qualifications)
    contracts <- lapply(
        record$workloads,
        metabPublicationValidateManifestWorkload
    )
    ids <- vapply(record$workloads, `[[`, character(1), "workload_id")
    classes <- vapply(record$workloads, `[[`, character(1), "workload_class")
    routes <- vapply(record$workloads, `[[`, character(1), "route")
    profiles <- vapply(record$workloads, `[[`, character(1), "profile_id")
    qualification <- stats::setNames(
        record$pilot_qualifications,
        vapply(record$pilot_qualifications, `[[`, character(1), "route")
    )
    heavy_valid <- all(vapply(seq_along(contracts), function(index) {
        if (!identical(classes[[index]], "operational_heavy")) return(TRUE)
        selected <- qualification[[routes[[index]]]]
        identical(
            contracts[[index]]$dimensions$aggregate_feature_count,
            selected$aggregate_feature_count
        ) && identical(
            contracts[[index]]$dimensions$sample_count,
            selected$sample_count
        )
    }, logical(1)))
    matrix_keys <- paste(routes, profiles, classes, sep = "::")
    valid <- identical(
        record$schema,
        "multischolar.omics_publication_metabolomics_manifest"
    ) && identical(record$schema_version, .METAB_PUBLICATION_VERSION) &&
        identical(record$owner_ticket_id, .METAB_PUBLICATION_OWNER) &&
        identical(record$status, "frozen_nonpromotional") &&
        length(record$workloads) == 32L && !anyDuplicated(ids) &&
        !anyDuplicated(matrix_keys) && all(vapply(
            metabPublicationClasses(),
            function(workload_class) sum(classes == workload_class) == 8L,
            logical(1)
        )) && heavy_valid &&
        identical(record$source_claim_scope, "project_specific_nonpromotional") &&
        isTRUE(record$candidate_access_allowed) &&
        !isTRUE(record$promotion_authority) &&
        !isTRUE(record$publication_authority)
    if (!valid) metabPublicationAbort("metabolomics final manifest differs")
    invisible(record)
}

metabPublicationValidateFinalHandoff <- function(record) {
    metabPublicationValidateHandoff(record)
    expected <- c(
        "manifest", "governance_successor", "sources", "splits",
        "support_boundaries", "mapping_authorities", "msdial_bundles",
        "exclusions", "negative", "pilot_qualifications"
    )
    if (!identical(names(record$producer_bindings), expected)) {
        metabPublicationAbort("metabolomics final handoff bindings differ")
    }
    invisible(record)
}
