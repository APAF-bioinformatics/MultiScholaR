proteomicsPublicationParameterValues <- function(authority) {
    proteomicsPublicationValidateParameters(authority)
    values <- lapply(authority$parameters, `[[`, "value")
    names(values) <- vapply(
        authority$parameters,
        `[[`,
        character(1),
        "parameter_id"
    )
    values
}

proteomicsPublicationParameter <- function(parameters, id) {
    value <- parameters[[id]]
    if (is.null(value)) {
        proteomicsPublicationAbort(paste("model parameter is missing:", id))
    }
    value
}

proteomicsPublicationSetSeed <- function(seed) {
    do.call(RNGkind, list("L'Ecuyer-CMRG", "Inversion", "Rejection"))
    set.seed(as.integer(seed))
    invisible(seed)
}

proteomicsPublicationBlockSeed <- function(seed, block_index) {
    modulus <- .Machine$integer.max - 1L
    value <- (as.double(seed) + as.double(block_index) * 104729) %% modulus
    as.integer(value + 1)
}

proteomicsPublicationSoftplus <- function(value) {
    log1p(exp(-abs(value))) + pmax(value, 0)
}

proteomicsPublicationSampleDesign <- function(contract, parameters) {
    sample_count <- as.integer(contract$dimensions$sample_count)
    treatment_fraction <- proteomicsPublicationParameter(
        parameters,
        "treatment_fraction"
    )
    treatment_count <- as.integer(round(sample_count * treatment_fraction))
    treatment_count <- max(2L, min(sample_count - 2L, treatment_count))
    control_count <- sample_count - treatment_count
    sample_ids <- c(
        sprintf("CTRL_%03d", seq_len(control_count)),
        sprintf("TREAT_%03d", seq_len(treatment_count))
    )
    batch_count <- as.integer(proteomicsPublicationParameter(
        parameters,
        "batch_count"
    ))
    technical_size <- as.integer(proteomicsPublicationParameter(
        parameters,
        "technical_replicate_size"
    ))
    data.frame(
        sample_id = sample_ids,
        group = c(rep("control", control_count), rep("treatment", treatment_count)),
        treated = c(rep(0, control_count), rep(1, treatment_count)),
        batch = sprintf("B%02d", (seq_len(sample_count) - 1L) %% batch_count + 1L),
        technical_group = sprintf(
            "TR%03d",
            (seq_len(sample_count) - 1L) %/% technical_size + 1L
        ),
        stringsAsFactors = FALSE
    )
}

proteomicsPublicationCorrelation <- function(design, parameters) {
    sample_count <- nrow(design)
    sample_rho <- proteomicsPublicationParameter(parameters, "sample_correlation")
    technical_rho <- proteomicsPublicationParameter(
        parameters,
        "technical_replicate_correlation"
    )
    proteomicsPublicationRequireNumber(
        sample_rho,
        "sample correlation",
        minimum = 0,
        maximum = 0.95
    )
    proteomicsPublicationRequireNumber(
        technical_rho,
        "technical correlation",
        minimum = sample_rho,
        maximum = 0.99
    )
    correlation <- matrix(sample_rho, nrow = sample_count, ncol = sample_count)
    same_technical <- outer(
        design$technical_group,
        design$technical_group,
        `==`
    )
    correlation[same_technical] <- technical_rho
    diag(correlation) <- 1
    eigenvalues <- eigen(correlation, symmetric = TRUE, only.values = TRUE)$values
    if (min(eigenvalues) <= 1e-8) {
        proteomicsPublicationAbort("sample correlation matrix is not positive definite")
    }
    list(
        matrix = correlation,
        cholesky = chol(correlation),
        minimum_eigenvalue = min(eigenvalues)
    )
}

proteomicsPublicationEffectVector <- function(protein_count, parameters) {
    fraction <- proteomicsPublicationParameter(parameters, "effect_fraction")
    effect_log2 <- proteomicsPublicationParameter(parameters, "effect_log2")
    proteomicsPublicationRequireNumber(
        fraction,
        "effect fraction",
        minimum = 0.001,
        maximum = 0.24
    )
    effect_count <- max(1L, as.integer(floor(protein_count * fraction)))
    effects <- numeric(protein_count)
    effects[seq_len(effect_count)] <- effect_log2
    down <- seq.int(effect_count + 1L, 2L * effect_count)
    effects[down] <- -effect_log2
    list(
        values = effects,
        up = seq_len(effect_count),
        down = down,
        unaffected = seq.int(2L * effect_count + 1L, protein_count)
    )
}

proteomicsPublicationProteinTruth <- function(contract, parameters) {
    protein_count <- as.integer(contract$dimensions$protein_count)
    seed <- contract$rng$streams$hierarchy
    proteomicsPublicationSetSeed(seed)
    baseline <- stats::rnorm(
        protein_count,
        mean = proteomicsPublicationParameter(parameters, "base_log2_mean"),
        sd = proteomicsPublicationParameter(parameters, "base_log2_sd")
    )
    effects <- proteomicsPublicationEffectVector(protein_count, parameters)
    data.frame(
        protein_id = sprintf("SYNPROT%08d", seq_len(protein_count)),
        gene_id = sprintf("SYNGENE%08d", seq_len(protein_count)),
        baseline_log2 = baseline,
        effect_log2 = effects$values,
        effect_class = c(
            rep("up", length(effects$up)),
            rep("down", length(effects$down)),
            rep("unaffected", length(effects$unaffected))
        ),
        stringsAsFactors = FALSE
    )
}

proteomicsPublicationAdjustCounts <- function(counts, total, maximum) {
    difference <- as.integer(total - sum(counts))
    if (difference > 0L) {
        while (difference > 0L) {
            eligible <- which(counts < maximum)
            if (!length(eligible)) {
                proteomicsPublicationAbort("entity multiplicity maximum is too low")
            }
            take <- eligible[seq_len(min(length(eligible), difference))]
            counts[take] <- counts[take] + 1L
            difference <- as.integer(total - sum(counts))
        }
    } else if (difference < 0L) {
        while (difference < 0L) {
            eligible <- which(counts > 1L)
            if (!length(eligible)) {
                proteomicsPublicationAbort("entity multiplicity total is too low")
            }
            take <- eligible[seq_len(min(length(eligible), abs(difference)))]
            counts[take] <- counts[take] - 1L
            difference <- as.integer(total - sum(counts))
        }
    }
    counts
}

proteomicsPublicationEntityCounts <- function(
    protein_count,
    entity_count,
    parameters,
    seed
) {
    maximum <- as.integer(proteomicsPublicationParameter(
        parameters,
        "maximum_peptides_per_protein"
    ))
    probabilities <- unlist(proteomicsPublicationParameter(
        parameters,
        "peptide_multiplicity_probabilities"
    ), use.names = FALSE)
    if (length(probabilities) != maximum || any(probabilities < 0) ||
        abs(sum(probabilities) - 1) > 1e-12) {
        proteomicsPublicationAbort("peptide multiplicity probabilities differ")
    }
    if (entity_count < protein_count || entity_count > protein_count * maximum) {
        proteomicsPublicationAbort("entity count is outside multiplicity bounds")
    }
    proteomicsPublicationSetSeed(seed)
    counts <- sample.int(
        maximum,
        size = protein_count,
        replace = TRUE,
        prob = probabilities
    )
    proteomicsPublicationAdjustCounts(counts, entity_count, maximum)
}

proteomicsPublicationPeptideTruth <- function(contract, proteins, parameters) {
    peptide_count <- as.integer(contract$dimensions$peptide_count)
    if (peptide_count == 0L) return(NULL)
    counts <- proteomicsPublicationEntityCounts(
        nrow(proteins),
        peptide_count,
        parameters,
        contract$rng$streams$multiplicity
    )
    protein_index <- rep.int(seq_len(nrow(proteins)), counts)
    proteomicsPublicationSetSeed(contract$rng$streams$peptide_offsets)
    offsets <- stats::rnorm(
        peptide_count,
        sd = proteomicsPublicationParameter(parameters, "peptide_offset_sd")
    )
    data.frame(
        peptide_id = sprintf("SYNPEPTIDE%09dK", seq_len(peptide_count)),
        protein_index = protein_index,
        protein_id = proteins$protein_id[protein_index],
        peptide_offset_log2 = offsets,
        stringsAsFactors = FALSE
    )
}

proteomicsPublicationPrecursorTruth <- function(contract, peptides) {
    precursor_count <- as.integer(contract$dimensions$precursor_count)
    if (precursor_count == 0L) return(NULL)
    peptide_count <- nrow(peptides)
    if (precursor_count < peptide_count || precursor_count > peptide_count * 3L) {
        proteomicsPublicationAbort("precursor count is outside charge bounds")
    }
    counts <- rep.int(1L, peptide_count)
    extra <- precursor_count - peptide_count
    if (extra > 0L) {
        first <- seq_len(min(extra, peptide_count))
        counts[first] <- counts[first] + 1L
        extra <- precursor_count - sum(counts)
    }
    if (extra > 0L) counts[seq_len(extra)] <- counts[seq_len(extra)] + 1L
    peptide_index <- rep.int(seq_len(peptide_count), counts)
    charge <- sequence(counts) + 1L
    data.frame(
        precursor_id = sprintf(
            "%s_%d",
            peptides$peptide_id[peptide_index],
            charge
        ),
        peptide_index = peptide_index,
        protein_index = peptides$protein_index[peptide_index],
        charge = charge,
        stringsAsFactors = FALSE
    )
}

proteomicsPublicationBatchOffsets <- function(design, parameters, seed) {
    batches <- unique(design$batch)
    proteomicsPublicationSetSeed(seed)
    offsets <- stats::rnorm(
        length(batches),
        sd = proteomicsPublicationParameter(parameters, "batch_log2_sd")
    )
    stats::setNames(offsets, batches)[design$batch]
}

proteomicsPublicationModelPlan <- function(contract, parameter_authority) {
    proteomicsPublicationValidateWorkload(contract)
    parameters <- proteomicsPublicationParameterValues(parameter_authority)
    design <- proteomicsPublicationSampleDesign(contract, parameters)
    correlation <- proteomicsPublicationCorrelation(design, parameters)
    proteins <- proteomicsPublicationProteinTruth(contract, parameters)
    peptides <- proteomicsPublicationPeptideTruth(
        contract,
        proteins,
        parameters
    )
    precursors <- proteomicsPublicationPrecursorTruth(contract, peptides)
    batch_offsets <- proteomicsPublicationBatchOffsets(
        design,
        parameters,
        contract$rng$streams$batch
    )
    list(
        contract = contract,
        parameters = parameters,
        design = design,
        correlation = correlation,
        proteins = proteins,
        peptides = peptides,
        precursors = precursors,
        batch_offsets = as.numeric(batch_offsets)
    )
}

proteomicsPublicationEntityMap <- function(plan) {
    format <- plan$contract$capability$input_format
    if (identical(format, "diann")) {
        return(data.frame(
            entity_id = plan$precursors$precursor_id,
            peptide_id = plan$peptides$peptide_id[
                plan$precursors$peptide_index
            ],
            protein_index = plan$precursors$protein_index,
            offset_log2 = plan$peptides$peptide_offset_log2[
                plan$precursors$peptide_index
            ],
            charge = plan$precursors$charge,
            stringsAsFactors = FALSE
        ))
    }
    if (identical(format, "pd_tmt")) {
        return(data.frame(
            entity_id = plan$peptides$peptide_id,
            peptide_id = plan$peptides$peptide_id,
            protein_index = plan$peptides$protein_index,
            offset_log2 = plan$peptides$peptide_offset_log2,
            charge = 0L,
            stringsAsFactors = FALSE
        ))
    }
    data.frame(
        entity_id = plan$proteins$protein_id,
        peptide_id = NA_character_,
        protein_index = seq_len(nrow(plan$proteins)),
        offset_log2 = 0,
        charge = 0L,
        stringsAsFactors = FALSE
    )
}

proteomicsPublicationRandomMatrix <- function(
    rows,
    columns,
    seed,
    block_index,
    distribution = c("normal", "uniform")
) {
    distribution <- match.arg(distribution)
    proteomicsPublicationSetSeed(proteomicsPublicationBlockSeed(
        seed,
        block_index
    ))
    values <- if (identical(distribution, "normal")) {
        stats::rnorm(rows * columns)
    } else {
        stats::runif(rows * columns)
    }
    matrix(values, nrow = rows, ncol = columns)
}

proteomicsPublicationGenerateBlock <- function(plan, entity_index) {
    entities <- proteomicsPublicationEntityMap(plan)[entity_index, , drop = FALSE]
    proteins <- plan$proteins[entities$protein_index, , drop = FALSE]
    design <- plan$design
    block_index <- (min(entity_index) - 1L) %/%
        as.integer(plan$contract$generator$chunk_rows) + 1L
    mean_log2 <- proteins$baseline_log2 + entities$offset_log2
    treated_effect <- proteins$effect_log2 %o% design$treated
    location <- mean_log2 + treated_effect
    location <- sweep(location, 2L, plan$batch_offsets, "+")
    sigma <- proteomicsPublicationParameter(
        plan$parameters,
        "residual_sigma_floor"
    ) + proteomicsPublicationParameter(
        plan$parameters,
        "residual_sigma_slope"
    ) * proteomicsPublicationSoftplus(
        proteomicsPublicationParameter(
            plan$parameters,
            "heteroscedastic_reference_log2"
        ) - mean_log2
    )
    residual <- proteomicsPublicationRandomMatrix(
        nrow(entities),
        nrow(design),
        plan$contract$rng$streams$residual,
        block_index
    ) %*% plan$correlation$cholesky
    residual <- residual * sigma
    latent <- location + residual
    mar_linear <- proteomicsPublicationParameter(
        plan$parameters,
        "mar_logit_intercept"
    ) + proteomicsPublicationParameter(
        plan$parameters,
        "mar_treatment_log_odds"
    ) * rep(design$treated, each = nrow(entities))
    mar_probability <- matrix(
        stats::plogis(mar_linear),
        nrow = nrow(entities),
        ncol = nrow(design)
    )
    detection <- stats::plogis(
        (latent - proteomicsPublicationParameter(
            plan$parameters,
            "detection_midpoint_log2"
        )) / proteomicsPublicationParameter(
            plan$parameters,
            "detection_slope_log2"
        )
    )
    mar_uniform <- proteomicsPublicationRandomMatrix(
        nrow(entities),
        nrow(design),
        plan$contract$rng$streams$mar,
        block_index,
        "uniform"
    )
    mnar_uniform <- proteomicsPublicationRandomMatrix(
        nrow(entities),
        nrow(design),
        plan$contract$rng$streams$mnar,
        block_index,
        "uniform"
    )
    mar_missing <- mar_uniform < mar_probability
    mnar_missing <- !mar_missing & mnar_uniform > detection
    values <- round(2^latent, digits = 6L)
    values[mar_missing | mnar_missing] <- NA_real_
    list(
        entities = entities,
        values = values,
        latent_log2 = latent,
        mar_missing = mar_missing,
        mnar_missing = mnar_missing,
        detection_probability = detection,
        residual_sigma = sigma
    )
}
