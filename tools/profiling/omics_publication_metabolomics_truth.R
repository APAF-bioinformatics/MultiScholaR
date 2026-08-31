metabPublicationTruthState <- function(plan) {
    state <- new.env(parent = emptyenv())
    state$quantity_count <- 0
    state$observed_count <- 0
    state$quantity_sum <- 0
    state$quantity_min <- Inf
    state$quantity_max <- -Inf
    state$mechanisms <- c(MCAR = 0, MAR = 0, MNAR = 0, left_censored = 0)
    assays <- unique(plan$features$assay)
    state$assay_quantities <- stats::setNames(numeric(length(assays)), assays)
    state$assay_observed <- stats::setNames(numeric(length(assays)), assays)
    state$block_count <- 0L
    state
}

metabPublicationObserveTruth <- function(state, block, index) {
    observed <- is.finite(block$values)
    values <- block$values[observed]
    state$quantity_count <- state$quantity_count + length(block$values)
    state$observed_count <- state$observed_count + length(values)
    state$quantity_sum <- state$quantity_sum + sum(values)
    if (length(values)) {
        state$quantity_min <- min(state$quantity_min, min(values))
        state$quantity_max <- max(state$quantity_max, max(values))
    }
    state$mechanisms <- state$mechanisms + c(
        MCAR = sum(block$mcar_missing),
        MAR = sum(block$mar_missing),
        MNAR = sum(block$mnar_missing),
        left_censored = sum(block$censored_missing)
    )
    for (assay in unique(block$features$assay)) {
        selected <- block$features$assay == assay
        state$assay_quantities[[assay]] <- state$assay_quantities[[assay]] +
            sum(selected) * ncol(block$values)
        state$assay_observed[[assay]] <- state$assay_observed[[assay]] +
            sum(observed[selected, , drop = FALSE])
    }
    state$block_count <- as.integer(index)
    invisible(state)
}

metabPublicationEffectTruth <- function(plan) {
    classes <- plan$features$effect_class
    list(
        effect_log2 = metabPublicationParameter(plan$parameters, "effect_log2"),
        minimum_sign_agreement = metabPublicationParameter(
            plan$parameters,
            "minimum_effect_sign_agreement"
        ),
        median_margin_fraction = metabPublicationParameter(
            plan$parameters,
            "effect_median_margin_fraction"
        ),
        up_count = as.integer(sum(classes == "up")),
        down_count = as.integer(sum(classes == "down")),
        unaffected_count = as.integer(sum(classes == "unaffected")),
        up_feature_ids = as.list(which(classes == "up")),
        down_feature_ids = as.list(which(classes == "down")),
        assignment_sha256 = publicationObjectDigest(plan$features[c(
            "feature_id", "effect_class", "effect_log2"
        )])
    )
}

metabPublicationValidateFixtureTruth <- function(record, contract) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "owner_ticket_id", "truth_id",
        "workload_id", "workload_class", "capability", "assay_profile",
        "payload", "source_review_bindings", "dimensions", "expected_import",
        "effects", "design", "mapping", "oracle_method", "verified_stages",
        "limitations", "publication_authority"
    ), "Metabolomics fixture truth")
    lapply(record$source_review_bindings, function(binding) {
        metabPublicationValidateBinding(binding, "Metabolomics fixture review")
    })
    expected <- record$expected_import
    valid <- identical(
        record$schema,
        "multischolar.omics_publication_metabolomics_fixture_truth"
    ) && identical(record$schema_version, "1.0.0") &&
        identical(record$owner_ticket_id, "OMICS-ART-066") &&
        identical(record$workload_id, contract$workload_id) &&
        identical(record$workload_class, "fixture_correctness") &&
        identical(record$capability, contract$capability) &&
        identical(record$assay_profile, contract$assay_profile) &&
        identical(record$dimensions, contract$dimensions) &&
        identical(
            record$payload$payload_set_sha256,
            contract$expected_digests$payload_set_sha256
        ) && identical(
            expected$aggregate_feature_count,
            contract$dimensions$aggregate_feature_count
        ) && identical(expected$sample_count, contract$dimensions$sample_count) &&
        expected$quantity_count == contract$dimensions$quantity_count &&
        length(record$effects$up_feature_ids) > 0L &&
        length(record$effects$down_feature_ids) > 0L &&
        identical(
            record$oracle_method,
            "direct_raw_table_arithmetic_and_reviewed_e2e_authority"
        ) && identical(
            record$verified_stages,
            contract$capability$verified_stages
        ) && length(record$limitations) > 0L &&
        !isTRUE(record$publication_authority)
    metabPublicationRequireDigest(
        record$payload$payload_set_sha256,
        "Metabolomics fixture payload"
    )
    if (!valid) metabPublicationAbort("metabolomics fixture truth differs")
    invisible(record)
}

metabPublicationFinalizeTruth <- function(
    state,
    plan,
    contract,
    payload
) {
    if (state$observed_count == 0L || !is.finite(state$quantity_min) ||
        !is.finite(state$quantity_max)) {
        metabPublicationAbort("metabolomics truth state is incomplete")
    }
    missing_count <- sum(state$mechanisms)
    list(
        schema = "multischolar.omics_publication_metabolomics_truth",
        schema_version = "1.0.0",
        owner_ticket_id = "OMICS-ART-066",
        truth_id = paste0(contract$workload_id, ".truth"),
        workload_id = contract$workload_id,
        workload_class = contract$workload_class,
        capability = contract$capability,
        assay_profile = contract$assay_profile,
        contract_basis_sha256 = metabPublicationContractBasis(contract),
        parameter_authority = contract$parameter_authority,
        mapping_authority = contract$mapping_authority,
        support_boundary = contract$support_boundary,
        payload = payload,
        dimensions = contract$dimensions,
        hierarchy = list(
            feature_plan_sha256 = publicationObjectDigest(plan$features),
            correlated_group_count = as.integer(length(unique(
                plan$features$correlated_group_id
            ))),
            groups_cross_assays = FALSE,
            internal_standard_count = as.integer(sum(
                plan$features$internal_standard
            ))
        ),
        design = list(
            sample_count = as.integer(nrow(plan$design)),
            control_count = as.integer(sum(plan$design$group == "control")),
            treatment_count = as.integer(sum(plan$design$group == "treatment")),
            batch_count = as.integer(length(unique(plan$design$batch))),
            design_sha256 = publicationObjectDigest(plan$design)
        ),
        covariance = list(
            matrix_sha256 = publicationObjectDigest(plan$correlation$matrix),
            minimum_eigenvalue = plan$correlation$minimum_eigenvalue,
            positive_definite = TRUE
        ),
        effects = metabPublicationEffectTruth(plan),
        observations = list(
            quantity_count = as.numeric(state$quantity_count),
            observed_count = as.numeric(state$observed_count),
            missing_count = as.numeric(missing_count),
            quantity_sum = state$quantity_sum,
            quantity_min = state$quantity_min,
            quantity_max = state$quantity_max,
            assay_quantity_counts = as.list(state$assay_quantities),
            assay_observed_counts = as.list(state$assay_observed),
            block_count = state$block_count
        ),
        missingness = list(
            mechanism_counts = as.list(state$mechanisms),
            overlap_count = 0L,
            precedence = "MCAR_then_MAR_then_MNAR_then_left_censoring",
            masks_independently_seeded = TRUE
        ),
        expected_import = list(
            member_count = payload$member_count,
            aggregate_feature_count = contract$dimensions$aggregate_feature_count,
            assay_feature_counts = contract$dimensions$assay_feature_counts,
            sample_count = contract$dimensions$sample_count,
            quantity_count = contract$dimensions$quantity_count,
            quantity_na_count = as.numeric(missing_count),
            quantity_sum = state$quantity_sum,
            numerical_tolerance = max(1, abs(state$quantity_sum)) * 1e-12
        ),
        model = list(
            abundance = paste0(
                "assay_group_baseline_plus_feature_offset_effect_batch_",
                "and_correlated_residual"
            ),
            missingness = "independent_seeded_masks_with_declared_precedence",
            chemical_identity_claim = FALSE
        ),
        verified_stages = contract$capability$verified_stages,
        limitations = list(
            "Synthetic correlated groups are not chemical adduct annotations.",
            "Generated evidence cannot widen route support or promotion scope."
        ),
        promotion_authority = FALSE,
        publication_authority = FALSE
    )
}

metabPublicationWriteTruth <- function(record, path) {
    if (file.exists(path) || dir.exists(path)) {
        metabPublicationAbort("metabolomics truth output already exists")
    }
    publicationWriteJson(record, path)
    list(
        path = path,
        sha256 = publicationFileDigest(path),
        size_bytes = as.numeric(file.info(path)$size)
    )
}

metabPublicationValidateTruth <- function(record, contract) {
    checks <- c(
        schema = identical(
            record$schema,
            "multischolar.omics_publication_metabolomics_truth"
        ),
        schema_version = identical(record$schema_version, "1.0.0"),
        owner = identical(record$owner_ticket_id, "OMICS-ART-066"),
        workload_id = identical(record$workload_id, contract$workload_id),
        workload_class = identical(
            record$workload_class,
            contract$workload_class
        ),
        capability = identical(
            publicationObjectDigest(record$capability),
            publicationObjectDigest(contract$capability)
        ),
        assay_profile = identical(
            publicationObjectDigest(record$assay_profile),
            publicationObjectDigest(contract$assay_profile)
        ),
        contract_basis = identical(
            record$contract_basis_sha256,
            metabPublicationContractBasis(contract)
        ),
        dimensions = identical(
            publicationObjectDigest(record$dimensions),
            publicationObjectDigest(contract$dimensions)
        ),
        covariance = isTRUE(record$covariance$positive_definite) &&
            record$covariance$minimum_eigenvalue > 0,
        hierarchy = !isTRUE(record$hierarchy$groups_cross_assays),
        missingness = identical(record$missingness$overlap_count, 0L),
        verified_stages = identical(
            record$verified_stages,
            contract$capability$verified_stages
        ),
        promotion = !isTRUE(record$promotion_authority),
        publication = !isTRUE(record$publication_authority)
    )
    metabPublicationRequireDigest(
        record$payload$payload_set_sha256,
        "Metabolomics truth payload"
    )
    if (!all(checks)) {
        metabPublicationAbort(paste(
            "metabolomics truth differs:",
            paste(names(checks)[!checks], collapse = ", ")
        ))
    }
    invisible(record)
}
