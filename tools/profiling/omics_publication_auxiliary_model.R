auxPublicationParameterDefinitions <- function() {
    limitation <- paste0(
        "Declared synthetic design parameter; no empirical biological, ",
        "instrument, cohort, or service-behaviour claim."
    )
    define <- function(
        parameter_id,
        value,
        unit,
        minimum,
        maximum,
        surface_id,
        consumer
    ) {
        list(
            parameter_id = parameter_id,
            value = value,
            unit = unit,
            domain = list(
                type = if (is.integer(value)) {
                    "integer_interval"
                } else {
                    "numeric_interval"
                },
                minimum = minimum,
                maximum = maximum
            ),
            origin = "methodological_design",
            source_binding = NULL,
            surface_id = surface_id,
            consumer = consumer,
            allowed_claim_vocabulary = list(
                "declared_synthetic",
                "mechanistic_simulation"
            ),
            limitations = list(limitation)
        )
    }
    list(
        define(
            "phosphosite_protein_pool_size",
            10000L,
            "proteins",
            10L,
            1000000L,
            "phosphosite.api",
            "auxPublicationPhosphositeFeaturePlan"
        ),
        define(
            "phosphosite_peptide_length",
            15L,
            "amino_acids",
            7L,
            45L,
            "phosphosite.api",
            "auxPublicationPhosphositeFeaturePlan"
        ),
        define(
            "phosphosite_ambiguous_accession_fraction",
            0.08,
            "fraction_of_evidence_rows",
            0,
            0.4,
            "phosphosite.api",
            "auxPublicationPhosphositeFeaturePlan"
        ),
        define(
            "phosphosite_missing_accession_fraction",
            0.02,
            "fraction_of_evidence_rows",
            0,
            0.2,
            "phosphosite.api",
            "auxPublicationPhosphositeFeaturePlan"
        ),
        define(
            "phosphosite_contaminant_fraction",
            0.01,
            "fraction_of_evidence_rows",
            0,
            0.2,
            "phosphosite.api",
            "auxPublicationPhosphositeFeaturePlan"
        ),
        define(
            "phosphosite_reverse_fraction",
            0.01,
            "fraction_of_evidence_rows",
            0,
            0.2,
            "phosphosite.api",
            "auxPublicationPhosphositeFeaturePlan"
        ),
        define(
            "phosphosite_zero_abundance_fraction",
            0.04,
            "fraction_of_evidence_rows",
            0,
            0.5,
            "phosphosite.api",
            "auxPublicationPhosphositeFeaturePlan"
        ),
        define(
            "phosphosite_duplicate_fraction",
            0.05,
            "fraction_of_evidence_rows",
            0,
            0.5,
            "phosphosite.api",
            "auxPublicationPhosphositeFeaturePlan"
        ),
        define(
            "phosphosite_multisite_fraction",
            0.25,
            "fraction_of_evidence_rows",
            0,
            0.8,
            "phosphosite.api",
            "auxPublicationPhosphositeFeaturePlan"
        ),
        define(
            "phosphosite_high_confidence_fraction",
            0.7,
            "fraction_of_evidence_rows",
            0.1,
            1,
            "phosphosite.api",
            "auxPublicationPhosphositeFeaturePlan"
        ),
        define(
            "phosphosite_repeat_match_fraction",
            0.03,
            "fraction_of_proteins",
            0,
            0.3,
            "phosphosite.api",
            "auxPublicationPhosphositeProteinPlan"
        ),
        define(
            "multiomics_layer_count",
            3L,
            "layers",
            3L,
            12L,
            "multiomics.api",
            "auxPublicationMultiomicsPlan"
        ),
        define(
            "multiomics_sample_overlap_fraction",
            0.75,
            "fraction_of_samples",
            0.1,
            1,
            "multiomics.api",
            "auxPublicationMultiomicsPlan"
        ),
        define(
            "multiomics_missing_weight_fraction",
            0.02,
            "fraction_of_feature_weights",
            0,
            0.5,
            "multiomics.api",
            "auxPublicationMultiomicsPlan"
        ),
        define(
            "multiomics_duplicate_feature_fraction",
            0.03,
            "fraction_of_feature_ids",
            0,
            0.5,
            "multiomics.api",
            "auxPublicationMultiomicsPlan"
        ),
        define(
            "multiomics_tied_weight_fraction",
            0.08,
            "fraction_of_feature_weights",
            0,
            0.5,
            "multiomics.api",
            "auxPublicationMultiomicsPlan"
        ),
        define(
            "multiomics_response_pathway_count",
            250L,
            "pathways_per_service",
            5L,
            100000L,
            "multiomics.api",
            "auxPublicationResponsePlan"
        ),
        define(
            "multiomics_response_identifier_multiplicity",
            4L,
            "identifiers_per_pathway",
            1L,
            1000L,
            "multiomics.api",
            "auxPublicationResponsePlan"
        )
    )
}

auxPublicationValidateParameter <- function(parameter) {
    publicationRequireNames(parameter, c(
        "parameter_id", "value", "unit", "domain", "origin",
        "source_binding", "surface_id", "consumer",
        "allowed_claim_vocabulary", "limitations"
    ), "Auxiliary parameter")
    publicationRequireNames(
        parameter$domain,
        c("type", "minimum", "maximum"),
        "Auxiliary parameter domain"
    )
    expected <- Filter(function(value) {
        identical(value$parameter_id, parameter$parameter_id)
    }, auxPublicationParameterDefinitions())
    valid <- length(expected) == 1L && identical(
        publicationObjectDigest(parameter),
        publicationObjectDigest(expected[[1L]])
    ) &&
        parameter$value >= parameter$domain$minimum &&
        parameter$value <= parameter$domain$maximum &&
        parameter$surface_id %in% c("phosphosite.api", "multiomics.api") &&
        publicationScalarString(parameter$consumer) &&
        identical(parameter$origin, "methodological_design") &&
        is.null(parameter$source_binding)
    if (!valid) auxPublicationAbort("auxiliary parameter differs")
    invisible(parameter)
}

auxPublicationBuildParameterAuthority <- function() {
    parameters <- auxPublicationParameterDefinitions()
    list(
        schema = "multischolar.omics_publication_auxiliary_parameters",
        schema_version = .AUX_PUBLICATION_VERSION,
        parameters_id = paste0(
            "multischolar.omics_publication_auxiliary_parameters.2026-08-28.v1"
        ),
        owner_ticket_id = .AUX_PUBLICATION_OWNER,
        status = "frozen_pre_candidate",
        parameters = parameters,
        consumer_registry = lapply(parameters, function(parameter) {
            list(
                parameter_id = parameter$parameter_id,
                consumer = parameter$consumer
            )
        }),
        generated_counts_as_real = FALSE,
        scientific_authority = FALSE,
        publication_authority = FALSE
    )
}

auxPublicationValidateParameters <- function(record) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "parameters_id", "owner_ticket_id",
        "status", "parameters", "consumer_registry",
        "generated_counts_as_real", "scientific_authority",
        "publication_authority"
    ), "Auxiliary parameter authority")
    lapply(record$parameters, auxPublicationValidateParameter)
    ids <- vapply(record$parameters, `[[`, character(1), "parameter_id")
    registry_ids <- vapply(
        record$consumer_registry,
        `[[`,
        character(1),
        "parameter_id"
    )
    registry_valid <- all(vapply(record$consumer_registry, function(value) {
        publicationRequireNames(
            value,
            c("parameter_id", "consumer"),
            "Auxiliary parameter consumer"
        )
        expected <- Filter(function(parameter) {
            identical(parameter$parameter_id, value$parameter_id)
        }, record$parameters)
        length(expected) == 1L &&
            identical(value$consumer, expected[[1L]]$consumer)
    }, logical(1)))
    expected_ids <- vapply(
        auxPublicationParameterDefinitions(),
        `[[`,
        character(1),
        "parameter_id"
    )
    valid <- identical(
        record$schema,
        "multischolar.omics_publication_auxiliary_parameters"
    ) && identical(record$schema_version, .AUX_PUBLICATION_VERSION) &&
        identical(
            record$parameters_id,
            "multischolar.omics_publication_auxiliary_parameters.2026-08-28.v1"
        ) &&
        identical(record$owner_ticket_id, .AUX_PUBLICATION_OWNER) &&
        identical(record$status, "frozen_pre_candidate") &&
        identical(ids, expected_ids) && !anyDuplicated(ids) &&
        identical(registry_ids, ids) && registry_valid &&
        !isTRUE(record$generated_counts_as_real) &&
        !isTRUE(record$scientific_authority) &&
        !isTRUE(record$publication_authority)
    if (!valid) auxPublicationAbort("auxiliary parameter authority differs")
    invisible(record)
}

auxPublicationParameter <- function(parameters, parameter_id) {
    matches <- Filter(function(parameter) {
        identical(parameter$parameter_id, parameter_id)
    }, parameters$parameters)
    if (length(matches) != 1L) {
        auxPublicationAbort("auxiliary parameter is missing or duplicated")
    }
    matches[[1L]]$value
}

auxPublicationFractionMask <- function(index, fraction, offset = 0L) {
    bucket <- (as.double(index) * 7919 + as.double(offset)) %% 10000
    bucket < round(fraction * 10000)
}

auxPublicationPhosphositeProteinPlan <- function(parameters, feature_count) {
    pool_limit <- as.integer(auxPublicationParameter(
        parameters,
        "phosphosite_protein_pool_size"
    ))
    protein_count <- min(pool_limit, max(10L, as.integer(ceiling(
        feature_count / 10
    ))))
    index <- seq_len(protein_count)
    peptide <- "AAAAAASTYAAAAAA"
    repeated <- auxPublicationFractionMask(
        index,
        auxPublicationParameter(
            parameters,
            "phosphosite_repeat_match_fraction"
        ),
        1103L
    )
    sequence <- rep.int(
        paste0("MMMMMMM", peptide, "GGGGGGG"),
        protein_count
    )
    sequence[repeated] <- paste0(
        "MMMMMMM",
        peptide,
        "GGGGGGG",
        peptide,
        "KKKKKKK"
    )
    accession <- paste0("SYN", sprintf("%09d", index))
    data.frame(
        db = "sp",
        uniprot_acc = accession,
        uniprot_id = paste0(accession, "_SYNTHETIC"),
        species = "Synthetic organism",
        tax_id = "0",
        gene_name = paste0("SYNGENE", sprintf("%09d", index)),
        protein_evidence = 1L,
        sequence_version = 1L,
        is_isoform = "Canonical",
        isoform_num = 0L,
        cleaned_acc = accession,
        status = "Reviewed",
        annotation_score = 5L,
        seq = sequence,
        seq_length = nchar(sequence, type = "chars"),
        repeated_peptide = repeated,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
}

auxPublicationPhosphositeFeaturePlan <- function(
    parameters,
    feature_count,
    sample_count,
    seed
) {
    index <- seq_len(feature_count)
    duplicate_fraction <- auxPublicationParameter(
        parameters,
        "phosphosite_duplicate_fraction"
    )
    duplicate <- auxPublicationFractionMask(index, duplicate_fraction, 2371L)
    base_index <- index
    base_index[duplicate & index > 1L] <- base_index[duplicate & index > 1L] - 1L
    proteins <- auxPublicationPhosphositeProteinPlan(parameters, feature_count)
    protein_index <- ((base_index - 1L) %% nrow(proteins)) + 1L
    category <- (as.double(base_index) * 3571 + as.double(seed)) %% 10000
    missing_limit <- 10000 * auxPublicationParameter(
        parameters,
        "phosphosite_missing_accession_fraction"
    )
    contaminant_limit <- missing_limit + 10000 * auxPublicationParameter(
        parameters,
        "phosphosite_contaminant_fraction"
    )
    reverse_limit <- contaminant_limit + 10000 * auxPublicationParameter(
        parameters,
        "phosphosite_reverse_fraction"
    )
    ambiguous_limit <- reverse_limit + 10000 * auxPublicationParameter(
        parameters,
        "phosphosite_ambiguous_accession_fraction"
    )
    category_id <- ifelse(
        category < missing_limit,
        "missing",
        ifelse(
            category < contaminant_limit,
            "contaminant",
            ifelse(
                category < reverse_limit,
                "reverse",
                ifelse(category < ambiguous_limit, "ambiguous", "single")
            )
        )
    )
    accession <- proteins$uniprot_acc[protein_index]
    second_index <- (protein_index %% nrow(proteins)) + 1L
    accession[category_id == "ambiguous"] <- paste(
        accession[category_id == "ambiguous"],
        proteins$uniprot_acc[second_index[category_id == "ambiguous"]],
        sep = ";"
    )
    accession[category_id == "missing"] <- NA_character_
    accession[category_id == "contaminant"] <- paste0(
        "CON__",
        accession[category_id == "contaminant"]
    )
    accession[category_id == "reverse"] <- paste0(
        "REV__",
        accession[category_id == "reverse"]
    )
    multisite <- auxPublicationFractionMask(
        base_index,
        auxPublicationParameter(parameters, "phosphosite_multisite_fraction"),
        3253L
    )
    high_confidence <- auxPublicationFractionMask(
        base_index,
        auxPublicationParameter(
            parameters,
            "phosphosite_high_confidence_fraction"
        ),
        4519L
    )
    primary_probability <- ifelse(high_confidence, 0.9, 0.6)
    phosphopeptide <- ifelse(
        multisite,
        sprintf("AAAAAAS(%.3f)T(0.800)YAAAAAA", primary_probability),
        sprintf("AAAAAAS(%.3f)TYAAAAAA", primary_probability)
    )
    zero_abundance <- auxPublicationFractionMask(
        base_index,
        auxPublicationParameter(
            parameters,
            "phosphosite_zero_abundance_fraction"
        ),
        5741L
    )
    sample_index <- seq_len(sample_count)
    values <- outer(
        as.double(base_index),
        as.double(sample_index),
        function(feature, sample) 1000 + feature * 0.01 + sample * 10
    )
    values[zero_abundance, ] <- 0
    colnames(values) <- paste0(
        "reporter intensity corrected_",
        sample_index
    )
    evidence <- data.frame(
        leading_proteins = accession,
        phospho_sty_probabilities = phosphopeptide,
        phospho_sty = ifelse(multisite, 2L, 1L),
        sequence = "AAAAAASTYAAAAAA",
        experiment = paste0("experiment_", ((base_index - 1L) %% 4L) + 1L),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    evidence <- cbind(
        evidence,
        as.data.frame(values, check.names = FALSE)
    )
    list(
        evidence = evidence,
        proteins = proteins,
        row_plan = data.frame(
            source_index = index,
            base_index = base_index,
            protein_index = protein_index,
            second_protein_index = second_index,
            category_id = category_id,
            duplicate = duplicate,
            multisite = multisite,
            high_confidence = high_confidence,
            zero_abundance = zero_abundance,
            stringsAsFactors = FALSE
        )
    )
}

auxPublicationLayerSampleRegistry <- function(
    layer_names,
    sample_count,
    overlap_fraction
) {
    shared_count <- max(1L, as.integer(floor(sample_count * overlap_fraction)))
    shared <- paste0("SYN_SAMPLE_SHARED_", sprintf("%05d", seq_len(shared_count)))
    stats::setNames(lapply(seq_along(layer_names), function(layer_index) {
        unique_count <- sample_count - shared_count
        unique_ids <- if (unique_count) {
            paste0(
                "SYN_SAMPLE_LAYER_",
                layer_index,
                "_",
                sprintf("%05d", seq_len(unique_count))
            )
        } else {
            character()
        }
        c(shared, unique_ids)
    }), layer_names)
}

auxPublicationMultiomicsPlan <- function(
    parameters,
    feature_count,
    sample_count,
    seed
) {
    layer_count <- as.integer(auxPublicationParameter(
        parameters,
        "multiomics_layer_count"
    ))
    layer_names <- c("proteome", "metabolome_lc", "metabolome_gc")
    if (layer_count != length(layer_names)) {
        auxPublicationAbort("multiomics layer count differs from API contract")
    }
    layer_index <- ((seq_len(feature_count) - 1L) %% layer_count) + 1L
    feature_index <- ave(seq_len(feature_count), layer_index, FUN = seq_along)
    missing <- auxPublicationFractionMask(
        seq_len(feature_count),
        auxPublicationParameter(
            parameters,
            "multiomics_missing_weight_fraction"
        ),
        6197L + seed
    )
    duplicated <- auxPublicationFractionMask(
        seq_len(feature_count),
        auxPublicationParameter(
            parameters,
            "multiomics_duplicate_feature_fraction"
        ),
        7151L + seed
    )
    tied <- auxPublicationFractionMask(
        seq_len(feature_count),
        auxPublicationParameter(parameters, "multiomics_tied_weight_fraction"),
        8089L + seed
    )
    feature_id <- paste0(
        "SYN_FEATURE_",
        sprintf("%09d", feature_index),
        "_",
        layer_names[layer_index]
    )
    feature_id[duplicated & feature_index > 1L] <- paste0(
        "SYN_FEATURE_",
        sprintf("%09d", feature_index[duplicated & feature_index > 1L] - 1L),
        "_",
        layer_names[layer_index[duplicated & feature_index > 1L]]
    )
    weight <- sin(as.double(seq_len(feature_count)) * 0.017 + seed / 1000)
    weight[tied] <- 0.5
    weight[missing] <- NA_real_
    weights <- lapply(seq_along(layer_names), function(layer) {
        selected <- layer_index == layer
        matrix(
            weight[selected],
            ncol = 1L,
            dimnames = list(feature_id[selected], "Factor1")
        )
    })
    names(weights) <- layer_names
    samples <- auxPublicationLayerSampleRegistry(
        layer_names,
        sample_count,
        auxPublicationParameter(
            parameters,
            "multiomics_sample_overlap_fraction"
        )
    )
    metabolite <- layer_index %in% c(2L, 3L)
    metabolite_index <- seq_len(sum(metabolite))
    metabolite_id <- feature_id[metabolite]
    view <- layer_names[layer_index[metabolite]]
    metabolite_base_id <- mapply(
        function(identifier, view_id) {
            sub(paste0("_", view_id, "$"), "", identifier)
        },
        metabolite_id,
        view,
        USE.NAMES = FALSE
    )
    chebi <- paste0("CHEBI:", 1000000L + metabolite_index)
    assay <- function(view_id) {
        selected <- view == view_id
        data.frame(
            metabolite = metabolite_base_id[selected],
            database_identifier = chebi[selected],
            metabolite_identification = metabolite_base_id[selected],
            SAMPLE_MEASURED_001 = 1000 + metabolite_index[selected],
            stringsAsFactors = FALSE,
            check.names = FALSE
        )
    }
    list(
        model_double = list(weights = weights),
        weights = data.frame(
            feature = feature_id,
            view = layer_names[layer_index],
            factor = "Factor1",
            value = weight,
            stringsAsFactors = FALSE
        ),
        assay_data = list(
            metabolome_lc = assay("metabolome_lc"),
            metabolome_gc = assay("metabolome_gc")
        ),
        mapping_table = data.frame(
            KEGG = paste0("C", sprintf("%05d", metabolite_index %% 100000L)),
            ChEBI = as.character(1000000L + metabolite_index),
            stringsAsFactors = FALSE
        ),
        sample_registry = samples,
        measured_sample_ids = "SAMPLE_MEASURED_001",
        feature_plan = data.frame(
            feature_id = feature_id,
            layer = layer_names[layer_index],
            missing = missing,
            duplicated = duplicated,
            tied = tied,
            stringsAsFactors = FALSE
        )
    )
}

auxPublicationResponsePlan <- function(
    parameters,
    responses,
    service_ids,
    row_count
) {
    multiplicity <- as.integer(auxPublicationParameter(
        parameters,
        "multiomics_response_identifier_multiplicity"
    ))
    pathway_count <- as.integer(auxPublicationParameter(
        parameters,
        "multiomics_response_pathway_count"
    ))
    stats::setNames(lapply(service_ids, function(service_id) {
        definition <- auxPublicationResponseDefinition(responses, service_id)
        auxPublicationExpandResponse(
            definition,
            row_count,
            multiplicity,
            pathway_count
        )
    }), service_ids)
}
