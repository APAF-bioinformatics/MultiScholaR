lipidPublicationParameter <- function(parameters, id) {
    value <- parameters[[id]]
    if (is.null(value)) lipidPublicationAbort(paste("parameter is missing:", id))
    value
}

lipidPublicationSetSeed <- function(seed) {
    do.call(RNGkind, list("L'Ecuyer-CMRG", "Inversion", "Rejection"))
    set.seed(as.integer(seed))
    invisible(seed)
}

lipidPublicationBlockSeed <- function(seed, block_index) {
    modulus <- .Machine$integer.max - 1L
    as.integer((as.double(seed) + block_index * 104729) %% modulus + 1)
}

lipidPublicationSoftplus <- function(value) {
    log1p(exp(-abs(value))) + pmax(value, 0)
}

lipidPublicationDesign <- function(sample_count, parameters) {
    treatment_fraction <- lipidPublicationParameter(
        parameters,
        "treatment_fraction"
    )
    treatment_count <- as.integer(round(sample_count * treatment_fraction))
    treatment_count <- max(3L, min(sample_count - 3L, treatment_count))
    control_count <- sample_count - treatment_count
    batch_count <- as.integer(lipidPublicationParameter(
        parameters,
        "batch_count"
    ))
    technical_size <- as.integer(lipidPublicationParameter(
        parameters,
        "technical_replicate_size"
    ))
    data.frame(
        sample_id = sprintf("LIPID_S%03d", seq_len(sample_count)),
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

lipidPublicationCorrelation <- function(design, parameters) {
    background <- lipidPublicationParameter(parameters, "sample_correlation")
    technical <- lipidPublicationParameter(
        parameters,
        "technical_replicate_correlation"
    )
    if (background < 0 || background >= 0.95 || technical < background ||
        technical >= 0.99) {
        lipidPublicationAbort("sample correlation parameters differ")
    }
    count <- nrow(design)
    matrix <- matrix(background, nrow = count, ncol = count)
    same <- outer(design$technical_group, design$technical_group, `==`)
    matrix[same] <- technical
    diag(matrix) <- 1
    eigenvalues <- eigen(matrix, symmetric = TRUE, only.values = TRUE)$values
    if (min(eigenvalues) <= 1e-8) {
        lipidPublicationAbort("sample correlation matrix is not positive definite")
    }
    list(
        matrix = matrix,
        cholesky = chol(matrix),
        minimum_eigenvalue = min(eigenvalues),
        lipid_class_correlation = lipidPublicationParameter(
            parameters,
            "lipid_class_correlation"
        )
    )
}

lipidPublicationAssayMeans <- function(parameters) {
    c(
        LCMS_Pos = lipidPublicationParameter(parameters, "lcms_pos_mean_log2"),
        LCMS_Neg = lipidPublicationParameter(parameters, "lcms_neg_mean_log2"),
        GCMS = lipidPublicationParameter(parameters, "gcms_mean_log2")
    )
}

lipidPublicationCensorThresholds <- function(parameters) {
    c(
        LCMS_Pos = lipidPublicationParameter(parameters, "lcms_pos_censor_log2"),
        LCMS_Neg = lipidPublicationParameter(parameters, "lcms_neg_censor_log2"),
        GCMS = lipidPublicationParameter(parameters, "gcms_censor_log2")
    )
}

lipidPublicationEffectVector <- function(feature_count, parameters) {
    fraction <- lipidPublicationParameter(parameters, "effect_fraction")
    count <- max(1L, as.integer(floor(feature_count * fraction)))
    effects <- numeric(feature_count)
    effects[seq_len(count)] <- lipidPublicationParameter(
        parameters,
        "effect_log2"
    )
    down <- seq.int(count + 1L, 2L * count)
    effects[down] <- -lipidPublicationParameter(parameters, "effect_log2")
    list(values = effects, up = seq_len(count), down = down)
}

lipidPublicationClassDefinitions <- function() {
    list(
        LCMS_Pos = list(
            classes = c("SYN_PC", "SYN_PE", "SYN_TG"),
            adducts = c("[M+H]+", "[M+H]+", "[M+NH4]+"),
            ion_mode = "positive"
        ),
        LCMS_Neg = list(
            classes = c("SYN_PI", "SYN_PS", "SYN_FA"),
            adducts = c("[M-H]-", "[M-H]-", "[M-H]-"),
            ion_mode = "negative"
        ),
        GCMS = list(
            classes = c("SYN_FA", "SYN_FAME", "SYN_STEROL"),
            adducts = c("[M]+.", "[M]+.", "[M]+."),
            ion_mode = "EI"
        )
    )
}

lipidPublicationFeatureStructure <- function(assays, parameters) {
    definitions <- lipidPublicationClassDefinitions()
    classes <- character(length(assays))
    adducts <- character(length(assays))
    ion_modes <- character(length(assays))
    for (assay in unique(assays)) {
        selected <- which(assays == assay)
        definition <- definitions[[assay]]
        local <- seq_along(selected)
        index <- (local - 1L) %% length(definition$classes) + 1L
        classes[selected] <- definition$classes[index]
        adducts[selected] <- definition$adducts[index]
        ion_modes[selected] <- definition$ion_mode
    }
    class_key <- paste(assays, classes, sep = "::")
    class_id <- match(class_key, unique(class_key))
    family_size <- as.integer(lipidPublicationParameter(
        parameters,
        "family_size"
    ))
    family_key <- character(length(assays))
    for (key in unique(class_key)) {
        selected <- which(class_key == key)
        family_key[selected] <- paste(
            key,
            (seq_along(selected) - 1L) %/% family_size + 1L,
            sep = "::"
        )
    }
    family_id <- match(family_key, unique(family_key))
    isomer_pair_id <- rep(NA_character_, length(assays))
    fraction <- lipidPublicationParameter(parameters, "isomer_like_fraction")
    pair_index <- 0L
    for (key in unique(class_key)) {
        selected <- which(class_key == key)
        pair_count <- as.integer(floor(length(selected) * fraction / 2L))
        if (!pair_count) next
        paired <- selected[seq_len(pair_count * 2L)]
        for (pair in split(paired, rep(seq_len(pair_count), each = 2L))) {
            pair_index <- pair_index + 1L
            isomer_pair_id[pair] <- sprintf("SYNISO%08d", pair_index)
        }
    }
    data.frame(
        lipid_class = classes,
        adduct = adducts,
        ion_mode = ion_modes,
        class_id = class_id,
        composition_family_id = sprintf("SYNCOMP%08d", family_id),
        family_index = family_id,
        isomer_pair_id = isomer_pair_id,
        stringsAsFactors = FALSE
    )
}

lipidPublicationFeaturePlan <- function(
    assay_feature_counts,
    parameters,
    streams
) {
    counts <- unlist(assay_feature_counts, use.names = TRUE)
    assays <- rep(names(counts), counts)
    feature_count <- length(assays)
    structure <- lipidPublicationFeatureStructure(assays, parameters)
    lipidPublicationSetSeed(streams$hierarchy)
    means <- lipidPublicationAssayMeans(parameters)
    baseline_deviation <- stats::rnorm(
        max(structure$family_index),
        mean = 0,
        sd = lipidPublicationParameter(parameters, "base_log2_sd")
    )
    class_deviation <- stats::rnorm(
        max(structure$class_id),
        mean = 0,
        sd = lipidPublicationParameter(parameters, "class_offset_log2_sd")
    )
    family_deviation <- stats::rnorm(
        max(structure$family_index),
        mean = 0,
        sd = lipidPublicationParameter(parameters, "family_offset_log2_sd")
    )
    lipidPublicationSetSeed(streams$feature_offsets)
    offsets <- stats::rnorm(
        feature_count,
        sd = lipidPublicationParameter(parameters, "feature_offset_log2_sd")
    )
    effects <- lipidPublicationEffectVector(feature_count, parameters)
    stride <- max(1L, as.integer(round(
        1 / lipidPublicationParameter(parameters, "internal_standard_fraction")
    )))
    internal <- seq_len(feature_count) %% stride == 0L
    effects$values[internal] <- 0
    effect_class <- ifelse(
        effects$values > 0,
        "up",
        ifelse(effects$values < 0, "down", "unaffected")
    )
    cbind(data.frame(
        feature_id = sprintf("SYNLIPID%09d", seq_len(feature_count)),
        assay = assays,
        baseline_log2 = means[assays] +
            class_deviation[structure$class_id] +
            baseline_deviation[structure$family_index] +
            family_deviation[structure$family_index] + offsets,
        effect_log2 = effects$values,
        effect_class = effect_class,
        internal_standard = internal,
        stringsAsFactors = FALSE
    ), structure)
}

lipidPublicationBatchOffsets <- function(design, parameters, seed) {
    batches <- unique(design$batch)
    lipidPublicationSetSeed(seed)
    values <- stats::rnorm(
        length(batches),
        sd = lipidPublicationParameter(parameters, "batch_log2_sd")
    )
    as.numeric(stats::setNames(values, batches)[design$batch])
}

lipidPublicationModelPlan <- function(
    assay_feature_counts,
    sample_count,
    parameter_authority,
    streams,
    chunk_rows
) {
    parameters <- lipidPublicationParameterValues(parameter_authority)
    design <- lipidPublicationDesign(sample_count, parameters)
    correlation <- lipidPublicationCorrelation(design, parameters)
    features <- lipidPublicationFeaturePlan(
        assay_feature_counts,
        parameters,
        streams
    )
    lipidPublicationSetSeed(streams$class_residual)
    class_residuals <- matrix(
        stats::rnorm(max(features$class_id) * nrow(design)),
        nrow = max(features$class_id)
    ) %*% correlation$cholesky
    list(
        parameters = parameters,
        design = design,
        correlation = correlation,
        features = features,
        class_residuals = class_residuals,
        batch_offsets = lipidPublicationBatchOffsets(
            design,
            parameters,
            streams$batch
        ),
        streams = streams,
        chunk_rows = as.integer(chunk_rows)
    )
}

lipidPublicationRandomMatrix <- function(
    rows,
    columns,
    seed,
    block_index,
    distribution = c("normal", "uniform")
) {
    distribution <- match.arg(distribution)
    lipidPublicationSetSeed(lipidPublicationBlockSeed(seed, block_index))
    values <- if (identical(distribution, "normal")) {
        stats::rnorm(rows * columns)
    } else {
        stats::runif(rows * columns)
    }
    matrix(values, nrow = rows, ncol = columns)
}

lipidPublicationGenerateBlock <- function(plan, feature_index) {
    features <- plan$features[feature_index, , drop = FALSE]
    block_index <- (min(feature_index) - 1L) %/% plan$chunk_rows + 1L
    location <- features$baseline_log2 +
        features$effect_log2 %o% plan$design$treated
    location <- sweep(location, 2L, plan$batch_offsets, "+")
    sigma <- lipidPublicationParameter(
        plan$parameters,
        "residual_sigma_floor"
    ) + lipidPublicationParameter(
        plan$parameters,
        "residual_sigma_slope"
    ) * lipidPublicationSoftplus(
        lipidPublicationParameter(
            plan$parameters,
            "heteroscedastic_reference_log2"
        ) - features$baseline_log2
    )
    feature_residual <- lipidPublicationRandomMatrix(
        nrow(features),
        nrow(plan$design),
        plan$streams$residual,
        block_index
    ) %*% plan$correlation$cholesky
    rho <- plan$correlation$lipid_class_correlation
    if (rho < 0 || rho >= 0.95) {
        lipidPublicationAbort("lipid class correlation differs")
    }
    residual <- sqrt(rho) * plan$class_residuals[features$class_id, ] +
        sqrt(1 - rho) * feature_residual
    latent <- location + residual * sigma
    mcar <- lipidPublicationRandomMatrix(
        nrow(features),
        nrow(plan$design),
        plan$streams$mcar,
        block_index,
        "uniform"
    ) < lipidPublicationParameter(plan$parameters, "mcar_probability")
    mar_probability <- stats::plogis(
        lipidPublicationParameter(plan$parameters, "mar_logit_intercept") +
            lipidPublicationParameter(plan$parameters, "mar_batch_log_odds") *
                rep(plan$design$batch != "B01", each = nrow(features))
    )
    mar <- !mcar & lipidPublicationRandomMatrix(
        nrow(features),
        nrow(plan$design),
        plan$streams$mar,
        block_index,
        "uniform"
    ) < matrix(mar_probability, nrow = nrow(features))
    detection <- stats::plogis(
        (latent - lipidPublicationParameter(
            plan$parameters,
            "detection_midpoint_log2"
        )) / lipidPublicationParameter(
            plan$parameters,
            "detection_slope_log2"
        )
    )
    mnar <- !mcar & !mar & lipidPublicationRandomMatrix(
        nrow(features),
        nrow(plan$design),
        plan$streams$mnar,
        block_index,
        "uniform"
    ) > detection
    thresholds <- lipidPublicationCensorThresholds(plan$parameters)[features$assay]
    censored <- !mcar & !mar & !mnar & latent < thresholds
    values <- round(2^latent, digits = 6L)
    values[mcar | mar | mnar | censored] <- NA_real_
    list(
        features = features,
        values = values,
        latent_log2 = latent,
        residual_sigma = sigma,
        mcar_missing = mcar,
        mar_missing = mar,
        mnar_missing = mnar,
        censored_missing = censored,
        detection_probability = detection
    )
}
