metabPublicationParameter <- function(parameters, id) {
    value <- parameters[[id]]
    if (is.null(value)) metabPublicationAbort(paste("parameter is missing:", id))
    value
}

metabPublicationSetSeed <- function(seed) {
    do.call(RNGkind, list("L'Ecuyer-CMRG", "Inversion", "Rejection"))
    set.seed(as.integer(seed))
    invisible(seed)
}

metabPublicationBlockSeed <- function(seed, block_index) {
    modulus <- .Machine$integer.max - 1L
    as.integer((as.double(seed) + block_index * 104729) %% modulus + 1)
}

metabPublicationSoftplus <- function(value) {
    log1p(exp(-abs(value))) + pmax(value, 0)
}

metabPublicationDesign <- function(sample_count, parameters) {
    treatment_fraction <- metabPublicationParameter(
        parameters,
        "treatment_fraction"
    )
    treatment_count <- as.integer(round(sample_count * treatment_fraction))
    treatment_count <- max(3L, min(sample_count - 3L, treatment_count))
    control_count <- sample_count - treatment_count
    batch_count <- as.integer(metabPublicationParameter(
        parameters,
        "batch_count"
    ))
    technical_size <- as.integer(metabPublicationParameter(
        parameters,
        "technical_replicate_size"
    ))
    data.frame(
        sample_id = sprintf("METAB_S%03d", seq_len(sample_count)),
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

metabPublicationCorrelation <- function(design, parameters) {
    background <- metabPublicationParameter(parameters, "sample_correlation")
    technical <- metabPublicationParameter(
        parameters,
        "technical_replicate_correlation"
    )
    if (background < 0 || background >= 0.95 || technical < background ||
        technical >= 0.99) {
        metabPublicationAbort("sample correlation parameters differ")
    }
    count <- nrow(design)
    matrix <- matrix(background, nrow = count, ncol = count)
    same <- outer(design$technical_group, design$technical_group, `==`)
    matrix[same] <- technical
    diag(matrix) <- 1
    eigenvalues <- eigen(matrix, symmetric = TRUE, only.values = TRUE)$values
    if (min(eigenvalues) <= 1e-8) {
        metabPublicationAbort("sample correlation matrix is not positive definite")
    }
    list(
        matrix = matrix,
        cholesky = chol(matrix),
        minimum_eigenvalue = min(eigenvalues)
    )
}

metabPublicationAssayMeans <- function(parameters) {
    c(
        LCMS_Pos = metabPublicationParameter(parameters, "lcms_pos_mean_log2"),
        LCMS_Neg = metabPublicationParameter(parameters, "lcms_neg_mean_log2"),
        GCMS = metabPublicationParameter(parameters, "gcms_mean_log2")
    )
}

metabPublicationCensorThresholds <- function(parameters) {
    c(
        LCMS_Pos = metabPublicationParameter(parameters, "lcms_pos_censor_log2"),
        LCMS_Neg = metabPublicationParameter(parameters, "lcms_neg_censor_log2"),
        GCMS = metabPublicationParameter(parameters, "gcms_censor_log2")
    )
}

metabPublicationEffectVector <- function(feature_count, parameters) {
    fraction <- metabPublicationParameter(parameters, "effect_fraction")
    count <- max(1L, as.integer(floor(feature_count * fraction)))
    effects <- numeric(feature_count)
    effects[seq_len(count)] <- metabPublicationParameter(
        parameters,
        "effect_log2"
    )
    down <- seq.int(count + 1L, 2L * count)
    effects[down] <- -metabPublicationParameter(parameters, "effect_log2")
    list(values = effects, up = seq_len(count), down = down)
}

metabPublicationFeaturePlan <- function(
    assay_feature_counts,
    parameters,
    streams
) {
    counts <- unlist(assay_feature_counts, use.names = TRUE)
    assays <- rep(names(counts), counts)
    feature_count <- length(assays)
    group_size <- as.integer(metabPublicationParameter(
        parameters,
        "correlated_group_size"
    ))
    group_counts <- as.integer(ceiling(counts / group_size))
    offsets <- c(0L, head(cumsum(group_counts), -1L))
    group_index <- unlist(lapply(seq_along(counts), \(index) {
        count <- counts[[index]]
        offsets[[index]] + (seq_len(count) - 1L) %/% group_size + 1L
    }), use.names = FALSE)
    metabPublicationSetSeed(streams$hierarchy)
    means <- metabPublicationAssayMeans(parameters)
    baseline_deviation <- stats::rnorm(
        max(group_index),
        mean = 0,
        sd = metabPublicationParameter(parameters, "base_log2_sd")
    )
    group_deviation <- stats::rnorm(
        max(group_index),
        mean = 0,
        sd = metabPublicationParameter(parameters, "group_offset_log2_sd")
    )
    group_baselines <- means[rep(names(counts), group_counts)] +
        baseline_deviation + group_deviation
    metabPublicationSetSeed(streams$feature_offsets)
    offsets <- stats::rnorm(
        feature_count,
        sd = metabPublicationParameter(parameters, "feature_offset_log2_sd")
    )
    effects <- metabPublicationEffectVector(feature_count, parameters)
    stride <- max(1L, as.integer(round(
        1 / metabPublicationParameter(parameters, "internal_standard_fraction")
    )))
    data.frame(
        feature_id = sprintf("SYNMETAB%09d", seq_len(feature_count)),
        assay = assays,
        group_index = group_index,
        correlated_group_id = sprintf("SYNGRP%08d", group_index),
        baseline_log2 = group_baselines[group_index] + offsets,
        effect_log2 = effects$values,
        effect_class = ifelse(
            seq_len(feature_count) %in% effects$up,
            "up",
            ifelse(seq_len(feature_count) %in% effects$down, "down", "unaffected")
        ),
        internal_standard = seq_len(feature_count) %% stride == 0L,
        stringsAsFactors = FALSE
    )
}

metabPublicationBatchOffsets <- function(design, parameters, seed) {
    batches <- unique(design$batch)
    metabPublicationSetSeed(seed)
    values <- stats::rnorm(
        length(batches),
        sd = metabPublicationParameter(parameters, "batch_log2_sd")
    )
    as.numeric(stats::setNames(values, batches)[design$batch])
}

metabPublicationModelPlan <- function(
    assay_feature_counts,
    sample_count,
    parameter_authority,
    streams,
    chunk_rows
) {
    parameters <- metabPublicationParameterValues(parameter_authority)
    design <- metabPublicationDesign(sample_count, parameters)
    list(
        parameters = parameters,
        design = design,
        correlation = metabPublicationCorrelation(design, parameters),
        features = metabPublicationFeaturePlan(
            assay_feature_counts,
            parameters,
            streams
        ),
        batch_offsets = metabPublicationBatchOffsets(
            design,
            parameters,
            streams$batch
        ),
        streams = streams,
        chunk_rows = as.integer(chunk_rows)
    )
}

metabPublicationRandomMatrix <- function(
    rows,
    columns,
    seed,
    block_index,
    distribution = c("normal", "uniform")
) {
    distribution <- match.arg(distribution)
    metabPublicationSetSeed(metabPublicationBlockSeed(seed, block_index))
    values <- if (identical(distribution, "normal")) {
        stats::rnorm(rows * columns)
    } else {
        stats::runif(rows * columns)
    }
    matrix(values, nrow = rows, ncol = columns)
}

metabPublicationGenerateBlock <- function(plan, feature_index) {
    features <- plan$features[feature_index, , drop = FALSE]
    block_index <- (min(feature_index) - 1L) %/% plan$chunk_rows + 1L
    location <- features$baseline_log2 +
        features$effect_log2 %o% plan$design$treated
    location <- sweep(location, 2L, plan$batch_offsets, "+")
    sigma <- metabPublicationParameter(
        plan$parameters,
        "residual_sigma_floor"
    ) + metabPublicationParameter(
        plan$parameters,
        "residual_sigma_slope"
    ) * metabPublicationSoftplus(
        metabPublicationParameter(
            plan$parameters,
            "heteroscedastic_reference_log2"
        ) - features$baseline_log2
    )
    residual <- metabPublicationRandomMatrix(
        nrow(features),
        nrow(plan$design),
        plan$streams$residual,
        block_index
    ) %*% plan$correlation$cholesky
    latent <- location + residual * sigma
    mcar <- metabPublicationRandomMatrix(
        nrow(features),
        nrow(plan$design),
        plan$streams$mcar,
        block_index,
        "uniform"
    ) < metabPublicationParameter(plan$parameters, "mcar_probability")
    mar_probability <- stats::plogis(
        metabPublicationParameter(plan$parameters, "mar_logit_intercept") +
            metabPublicationParameter(plan$parameters, "mar_batch_log_odds") *
                rep(plan$design$batch != "B01", each = nrow(features))
    )
    mar <- !mcar & metabPublicationRandomMatrix(
        nrow(features),
        nrow(plan$design),
        plan$streams$mar,
        block_index,
        "uniform"
    ) < matrix(mar_probability, nrow = nrow(features))
    detection <- stats::plogis(
        (latent - metabPublicationParameter(
            plan$parameters,
            "detection_midpoint_log2"
        )) / metabPublicationParameter(
            plan$parameters,
            "detection_slope_log2"
        )
    )
    mnar <- !mcar & !mar & metabPublicationRandomMatrix(
        nrow(features),
        nrow(plan$design),
        plan$streams$mnar,
        block_index,
        "uniform"
    ) > detection
    thresholds <- metabPublicationCensorThresholds(plan$parameters)[features$assay]
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
