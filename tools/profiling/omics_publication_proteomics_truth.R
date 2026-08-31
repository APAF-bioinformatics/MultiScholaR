proteomicsPublicationTruthState <- function(plan) {
    state <- new.env(parent = emptyenv())
    state$quantity_count <- 0
    state$observed_count <- 0
    state$mar_missing_count <- 0
    state$mnar_missing_count <- 0
    state$quantity_sum <- 0
    state$quantity_compensation <- 0
    state$quantity_min <- Inf
    state$quantity_max <- -Inf
    state$latent_count <- 0
    state$latent_mean <- 0
    state$latent_m2 <- 0
    state$sigma_min <- Inf
    state$sigma_max <- -Inf
    state$sample_observed <- numeric(nrow(plan$design))
    state$block_count <- 0L
    state
}

proteomicsPublicationKahanAdd <- function(total, compensation, value) {
    adjusted <- value - compensation
    updated <- total + adjusted
    list(
        total = updated,
        compensation = (updated - total) - adjusted
    )
}

proteomicsPublicationCombineMoments <- function(
    count,
    mean,
    m2,
    block_count,
    block_mean,
    block_m2
) {
    if (block_count == 0L) return(list(count = count, mean = mean, m2 = m2))
    if (count == 0L) {
        return(list(count = block_count, mean = block_mean, m2 = block_m2))
    }
    difference <- block_mean - mean
    combined <- count + block_count
    list(
        count = combined,
        mean = mean + difference * block_count / combined,
        m2 = m2 + block_m2 + difference^2 * count * block_count / combined
    )
}

proteomicsPublicationObserveTruth <- function(state, block, index) {
    observed <- is.finite(block$values)
    values <- block$values[observed]
    state$quantity_count <- state$quantity_count + length(block$values)
    state$observed_count <- state$observed_count + length(values)
    state$mar_missing_count <- state$mar_missing_count + sum(block$mar_missing)
    state$mnar_missing_count <- state$mnar_missing_count + sum(block$mnar_missing)
    if (length(values)) {
        addition <- proteomicsPublicationKahanAdd(
            state$quantity_sum,
            state$quantity_compensation,
            sum(values)
        )
        state$quantity_sum <- addition$total
        state$quantity_compensation <- addition$compensation
        state$quantity_min <- min(state$quantity_min, min(values))
        state$quantity_max <- max(state$quantity_max, max(values))
    }
    latent <- as.numeric(block$latent_log2)
    block_mean <- mean(latent)
    moments <- proteomicsPublicationCombineMoments(
        state$latent_count,
        state$latent_mean,
        state$latent_m2,
        length(latent),
        block_mean,
        sum((latent - block_mean)^2)
    )
    state$latent_count <- moments$count
    state$latent_mean <- moments$mean
    state$latent_m2 <- moments$m2
    state$sigma_min <- min(state$sigma_min, block$residual_sigma)
    state$sigma_max <- max(state$sigma_max, block$residual_sigma)
    state$sample_observed <- state$sample_observed + colSums(observed)
    state$block_count <- as.integer(index)
    invisible(state)
}

proteomicsPublicationHierarchyTruth <- function(plan) {
    peptide_multiplicity <- if (is.null(plan$peptides)) {
        integer()
    } else {
        tabulate(
            plan$peptides$protein_index,
            nbins = nrow(plan$proteins)
        )
    }
    precursor_multiplicity <- if (is.null(plan$precursors)) {
        integer()
    } else {
        tabulate(
            plan$precursors$peptide_index,
            nbins = nrow(plan$peptides)
        )
    }
    list(
        protein_count = as.integer(nrow(plan$proteins)),
        peptide_count = as.integer(if (is.null(plan$peptides)) {
            0L
        } else {
            nrow(plan$peptides)
        }),
        precursor_count = as.integer(if (is.null(plan$precursors)) {
            0L
        } else {
            nrow(plan$precursors)
        }),
        protein_truth_sha256 = publicationObjectDigest(plan$proteins),
        peptide_truth_sha256 = publicationObjectDigest(plan$peptides),
        precursor_truth_sha256 = publicationObjectDigest(plan$precursors),
        peptide_multiplicity = list(
            minimum = if (length(peptide_multiplicity)) {
                min(peptide_multiplicity)
            } else {
                0L
            },
            maximum = if (length(peptide_multiplicity)) {
                max(peptide_multiplicity)
            } else {
                0L
            },
            sum = as.integer(sum(peptide_multiplicity)),
            sha256 = publicationObjectDigest(as.list(peptide_multiplicity))
        ),
        precursor_multiplicity = list(
            minimum = if (length(precursor_multiplicity)) {
                min(precursor_multiplicity)
            } else {
                0L
            },
            maximum = if (length(precursor_multiplicity)) {
                max(precursor_multiplicity)
            } else {
                0L
            },
            sum = as.integer(sum(precursor_multiplicity)),
            sha256 = publicationObjectDigest(as.list(precursor_multiplicity))
        )
    )
}

proteomicsPublicationEffectTruth <- function(plan) {
    up <- plan$proteins$protein_id[plan$proteins$effect_class == "up"]
    down <- plan$proteins$protein_id[plan$proteins$effect_class == "down"]
    unaffected <- plan$proteins$protein_id[
        plan$proteins$effect_class == "unaffected"
    ]
    list(
        effect_log2 = proteomicsPublicationParameter(
            plan$parameters,
            "effect_log2"
        ),
        minimum_sign_agreement = proteomicsPublicationParameter(
            plan$parameters,
            "minimum_effect_sign_agreement"
        ),
        median_margin_fraction = proteomicsPublicationParameter(
            plan$parameters,
            "effect_median_margin_fraction"
        ),
        up_count = as.integer(length(up)),
        down_count = as.integer(length(down)),
        unaffected_count = as.integer(length(unaffected)),
        up_ids_sha256 = publicationObjectDigest(as.list(up)),
        down_ids_sha256 = publicationObjectDigest(as.list(down)),
        unaffected_ids_sha256 = publicationObjectDigest(as.list(unaffected)),
        assignment_rule = "first_up_then_down_then_unaffected"
    )
}

proteomicsPublicationObservationTruth <- function(state) {
    if (state$observed_count == 0L || state$latent_count < 2L ||
        !is.finite(state$quantity_min) || !is.finite(state$quantity_max)) {
        proteomicsPublicationAbort("truth observation state is incomplete")
    }
    list(
        quantity_count = as.numeric(state$quantity_count),
        observed_count = as.numeric(state$observed_count),
        missing_count = as.numeric(
            state$mar_missing_count + state$mnar_missing_count
        ),
        quantity_sum = state$quantity_sum,
        quantity_min = state$quantity_min,
        quantity_max = state$quantity_max,
        latent_log2_mean = state$latent_mean,
        latent_log2_sd = sqrt(state$latent_m2 / (state$latent_count - 1)),
        residual_sigma_min = state$sigma_min,
        residual_sigma_max = state$sigma_max,
        sample_observed_counts = as.list(state$sample_observed),
        block_count = state$block_count
    )
}

proteomicsPublicationMissingnessTruth <- function(state, plan) {
    total <- state$quantity_count
    list(
        mechanisms = list(
            mar = list(
                definition = "logistic_on_observed_group_and_batch_covariates",
                missing_count = as.numeric(state$mar_missing_count),
                fraction = state$mar_missing_count / total
            ),
            intensity_dependent = list(
                definition = "logistic_detection_on_latent_log2_intensity",
                classification = "MNAR_left_censoring",
                missing_count = as.numeric(state$mnar_missing_count),
                fraction = state$mnar_missing_count / total
            )
        ),
        overlap_count = 0L,
        mask_precedence = "MAR_then_intensity_dependent_on_remaining_cells",
        parameter_authority_sha256 = plan$contract$parameter_authority$sha256
    )
}

proteomicsPublicationExpectedImport <- function(state, plan) {
    format <- plan$contract$capability$input_format
    data_type <- if (identical(format, "diann")) "peptide" else "protein"
    list(
        data_type = data_type,
        row_count = as.numeric(plan$contract$dimensions$quantity_count),
        protein_count = as.integer(plan$contract$dimensions$protein_count),
        peptide_count = as.integer(plan$contract$dimensions$peptide_count),
        sample_count = as.integer(plan$contract$dimensions$sample_count),
        quantity_count = as.numeric(plan$contract$dimensions$quantity_count),
        quantity_na_count = as.numeric(
            state$mar_missing_count + state$mnar_missing_count
        ),
        quantity_sum = state$quantity_sum,
        quantity_min = state$quantity_min,
        quantity_max = state$quantity_max,
        numerical_tolerance = max(1, abs(state$quantity_sum)) * 1e-12
    )
}

proteomicsPublicationFinalizeTruth <- function(state, plan, payload_binding) {
    list(
        schema = "multischolar.omics_publication_proteomics_truth",
        schema_version = "1.0.0",
        owner_ticket_id = "OMICS-ART-065",
        truth_id = paste0(plan$contract$workload_id, ".truth"),
        workload_id = plan$contract$workload_id,
        workload_class = plan$contract$workload_class,
        capability = plan$contract$capability,
        contract_basis_sha256 = proteomicsPublicationContractBasis(
            plan$contract
        ),
        parameter_authority = plan$contract$parameter_authority,
        source_authority = plan$contract$source_authority,
        split_authority = plan$contract$split_authority,
        payload = payload_binding,
        dimensions = plan$contract$dimensions,
        hierarchy = proteomicsPublicationHierarchyTruth(plan),
        design = list(
            sample_count = as.integer(nrow(plan$design)),
            control_count = as.integer(sum(plan$design$group == "control")),
            treatment_count = as.integer(sum(plan$design$group == "treatment")),
            batch_count = as.integer(length(unique(plan$design$batch))),
            technical_group_count = as.integer(length(unique(
                plan$design$technical_group
            ))),
            design_sha256 = publicationObjectDigest(plan$design)
        ),
        covariance = list(
            matrix_sha256 = publicationObjectDigest(plan$correlation$matrix),
            minimum_eigenvalue = plan$correlation$minimum_eigenvalue,
            positive_definite = TRUE
        ),
        effects = proteomicsPublicationEffectTruth(plan),
        observations = proteomicsPublicationObservationTruth(state),
        missingness = proteomicsPublicationMissingnessTruth(state, plan),
        expected_import = proteomicsPublicationExpectedImport(state, plan),
        model = list(
            abundance = paste0(
                "protein_baseline_plus_entity_offset_plus_effect_plus_batch_",
                "plus_correlated_residual"
            ),
            heteroscedasticity = "sigma_floor_plus_sigma_slope_times_softplus_reference_minus_mean",
            covariance = "validated_compound_sample_and_technical_block_correlation",
            value_scale = "round_2_power_latent_log2_to_6_decimals"
        ),
        limitations = list(
            "Declared synthetic truth; no instrument or cohort realism claim.",
            "Generated evidence cannot establish parser support or promotion."
        ),
        publication_authority = FALSE
    )
}

proteomicsPublicationWriteTruth <- function(record, path) {
    if (file.exists(path) || dir.exists(path)) {
        proteomicsPublicationAbort("truth output already exists")
    }
    publicationWriteJson(record, path)
    list(
        path = path,
        sha256 = publicationFileDigest(path),
        size_bytes = as.numeric(file.info(path)$size)
    )
}

proteomicsPublicationValidateTruth <- function(record, contract) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "owner_ticket_id", "truth_id",
        "workload_id", "workload_class", "capability",
        "contract_basis_sha256", "parameter_authority", "source_authority",
        "split_authority", "payload", "dimensions", "hierarchy", "design",
        "covariance", "effects", "observations", "missingness",
        "expected_import", "model", "limitations", "publication_authority"
    ), "Proteomics publication truth")
    expected <- record$expected_import
    valid <- identical(
        record$schema,
        "multischolar.omics_publication_proteomics_truth"
    ) && identical(record$schema_version, "1.0.0") &&
        identical(record$owner_ticket_id, "OMICS-ART-065") &&
        identical(record$workload_id, contract$workload_id) &&
        identical(record$workload_class, contract$workload_class) &&
        identical(record$capability, contract$capability) &&
        identical(
            record$contract_basis_sha256,
            proteomicsPublicationContractBasis(contract)
        ) && identical(record$parameter_authority, contract$parameter_authority) &&
        identical(record$source_authority, contract$source_authority) &&
        identical(record$split_authority, contract$split_authority) &&
        identical(record$dimensions, contract$dimensions) &&
        identical(record$hierarchy$protein_count, contract$dimensions$protein_count) &&
        identical(record$hierarchy$peptide_count, contract$dimensions$peptide_count) &&
        identical(
            record$hierarchy$precursor_count,
            contract$dimensions$precursor_count
        ) && isTRUE(record$covariance$positive_definite) &&
        record$covariance$minimum_eigenvalue > 0 &&
        expected$quantity_count == contract$dimensions$quantity_count &&
        expected$row_count == contract$dimensions$quantity_count &&
        identical(record$missingness$overlap_count, 0L) &&
        !isTRUE(record$publication_authority)
    for (field in c("sha256")) {
        proteomicsPublicationRequireDigest(record$payload[[field]], "truth payload")
    }
    if (!valid) proteomicsPublicationAbort("truth contract differs")
    invisible(record)
}

proteomicsPublicationValidateFixtureTruth <- function(record, contract) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "owner_ticket_id", "truth_id",
        "workload_id", "capability", "payload", "review_authority",
        "dimensions", "expected_import", "effects", "oracle_method",
        "limitations", "publication_authority"
    ), "Proteomics fixture truth")
    proteomicsPublicationValidateBinding(
        record$review_authority,
        "Fixture review authority"
    )
    valid <- identical(
        record$schema,
        "multischolar.omics_publication_proteomics_fixture_truth"
    ) && identical(record$schema_version, "1.0.0") &&
        identical(record$owner_ticket_id, "OMICS-ART-065") &&
        identical(record$workload_id, contract$workload_id) &&
        identical(record$capability, contract$capability) &&
        identical(record$payload$sha256, contract$expected_digests$payload_sha256) &&
        identical(record$dimensions, contract$dimensions) &&
        identical(record$expected_import$row_count, contract$dimensions$quantity_count) &&
        identical(
            record$expected_import$protein_count,
            contract$dimensions$protein_count
        ) && identical(
            record$expected_import$peptide_count,
            contract$dimensions$peptide_count
        ) && identical(
            record$expected_import$sample_count,
            contract$dimensions$sample_count
        ) && identical(
            record$oracle_method,
            "direct_raw_table_arithmetic_and_hand_reviewed_readme"
        ) && length(record$limitations) > 0L &&
        !isTRUE(record$publication_authority)
    proteomicsPublicationRequireDigest(record$payload$sha256, "Fixture payload")
    if (!valid) proteomicsPublicationAbort("fixture truth differs")
    invisible(record)
}
