.METAB_WORKLOAD_ADAPTER_ID <- "multischolar.metabolomics.custom.v1"
.METAB_WORKLOAD_ADAPTER_VERSION <- "1.0.0"

metabWorkloadAbort <- function(message) {
    stop(message, call. = FALSE)
}

metabWorkloadDigest <- function(value) {
    canonicalize <- function(candidate) {
        if (!is.list(candidate)) return(candidate)
        if (!is.null(names(candidate))) {
            candidate <- candidate[order(names(candidate), method = "radix")]
        }
        lapply(candidate, canonicalize)
    }
    digest::digest(
        jsonlite::toJSON(
            canonicalize(value),
            auto_unbox = TRUE,
            null = "null",
            na = "null",
            digits = 17
        ),
        algo = "sha256",
        serialize = FALSE
    )
}

metabWorkloadParameters <- function(contract) {
    parameters <- contract$adapter_parameters
    required <- c(
        "generator_version", "profile", "internal_standard_fraction",
        "effect_fraction", "effect_log2", "batch_log2",
        "heteroscedasticity", "bounded_query_features",
        "query_repetitions"
    )
    if (!is.list(parameters) || !setequal(names(parameters), required)) {
        metabWorkloadAbort(
            "metabolomics adapter parameters do not match the contract"
        )
    }
    if (!identical(parameters$generator_version, "1.0.0")) {
        metabWorkloadAbort("unsupported metabolomics generator version")
    }
    parameters
}

metabWorkloadAssays <- function(contract) {
    assays <- names(contract$assay_mix)
    supported <- c("LCMS_Pos", "LCMS_Neg", "GCMS")
    if (!length(assays) || !all(assays %in% supported) ||
        length(assays) != as.integer(contract$dimensions$assay_count)) {
        metabWorkloadAbort("metabolomics assay mix is unsupported")
    }
    assays
}

metabWorkloadSamples <- function(sample_count) {
    if (sample_count < 6L || sample_count %% 2L != 0L) {
        metabWorkloadAbort("metabolomics sample_count must be even and at least six")
    }
    group_size <- sample_count / 2L
    data.frame(
        sample_id = sprintf("METAB_S%02d", seq_len(sample_count)),
        group = rep(c("CTRL", "TREAT"), each = group_size),
        batch = rep(c("B1", "B2"), length.out = sample_count),
        stringsAsFactors = FALSE
    )
}

metabWorkloadAssayMeans <- function(contract, assays) {
    declared <- contract$distributions$assay_log2_means
    values <- unlist(declared, use.names = TRUE)
    if (!identical(sort(names(values)), sort(assays)) ||
        any(!is.finite(as.numeric(values)))) {
        metabWorkloadAbort("assay intensity means are incomplete")
    }
    as.numeric(values[assays])
}

metabWorkloadAssayAssignment <- function(feature_count, assays) {
    rep(assays, length.out = feature_count)
}

metabWorkloadEffects <- function(feature_count, parameters) {
    fraction <- as.numeric(parameters$effect_fraction)
    if (!is.finite(fraction) || fraction <= 0 || fraction >= 0.25) {
        metabWorkloadAbort("effect_fraction must be between zero and 0.25")
    }
    count <- max(1L, as.integer(floor(feature_count * fraction)))
    effects <- numeric(feature_count)
    effects[seq_len(count)] <- as.numeric(parameters$effect_log2)
    effects[seq.int(count + 1L, 2L * count)] <-
        -as.numeric(parameters$effect_log2)
    list(values = effects, count = count)
}

metabWorkloadIntensityMatrix <- function(
    contract,
    parameters,
    assays,
    assignment,
    design,
    effects
) {
    feature_count <- length(assignment)
    sample_count <- nrow(design)
    means <- stats::setNames(
        metabWorkloadAssayMeans(contract, assays),
        assays
    )
    base_sd <- as.numeric(contract$distributions$base_log2_sd)
    noise_sd <- as.numeric(contract$distributions$noise_log2_sd)
    hetero <- as.numeric(parameters$heteroscedasticity)
    baseline <- stats::rnorm(
        feature_count,
        mean = means[assignment],
        sd = base_sd
    )
    treated <- as.integer(design$group == "TREAT")
    batch <- ifelse(design$batch == "B2", as.numeric(parameters$batch_log2), 0)
    feature_noise <- noise_sd * (
        1 + hetero * abs(baseline - mean(baseline)) /
            max(stats::sd(baseline), .Machine$double.eps)
    )
    noise <- matrix(
        stats::rnorm(feature_count * sample_count),
        nrow = feature_count,
        ncol = sample_count
    ) * feature_noise
    log_values <- baseline + effects %o% treated
    log_values <- sweep(log_values, 2L, batch, "+") + noise
    round(pmax(2^log_values, 0), digits = 6L)
}

metabWorkloadMechanisms <- function(feature_count) {
    mechanisms <- c("MCAR", "MAR", "MNAR", "left_censored")
    rep(mechanisms, length.out = feature_count)
}

metabWorkloadApplyMissingness <- function(
    values,
    mechanisms,
    design,
    assignment,
    contract
) {
    counts <- stats::setNames(integer(4L), c(
        "MCAR", "MAR", "MNAR", "left_censored"
    ))
    for (row in seq_len(nrow(values))) {
        mechanism <- mechanisms[[row]]
        if (identical(mechanism, "MCAR")) {
            columns <- which((row + seq_len(ncol(values))) %% 19L == 0L)
        } else if (identical(mechanism, "MAR")) {
            columns <- which(design$batch == "B2" &
                (row + seq_len(ncol(values))) %% 7L == 0L)
        } else if (identical(mechanism, "MNAR")) {
            columns <- order(values[row, ], method = "radix")[[1L]]
        } else {
            thresholds <- unlist(
                contract$censoring$assay_intensity_thresholds,
                use.names = TRUE
            )
            threshold <- as.numeric(thresholds[[assignment[[row]]]])
            columns <- which(values[row, ] < threshold)
        }
        if (length(columns)) {
            values[row, columns] <- NA_real_
            counts[[mechanism]] <- counts[[mechanism]] + length(columns)
        }
    }
    list(values = values, counts = counts)
}

metabWorkloadPayload <- function(contract, parameters) {
    feature_count <- as.integer(contract$dimensions$feature_count)
    sample_count <- as.integer(contract$dimensions$sample_count)
    assays <- metabWorkloadAssays(contract)
    design <- metabWorkloadSamples(sample_count)
    assignment <- metabWorkloadAssayAssignment(feature_count, assays)
    effects <- metabWorkloadEffects(feature_count, parameters)
    values <- metabWorkloadIntensityMatrix(
        contract,
        parameters,
        assays,
        assignment,
        design,
        effects$values
    )
    mechanisms <- metabWorkloadMechanisms(feature_count)
    missing <- metabWorkloadApplyMissingness(
        values,
        mechanisms,
        design,
        assignment,
        contract
    )
    values <- missing$values
    values[[1L, 1L]] <- 0
    feature_ids <- sprintf("METAB_F%07d", seq_len(feature_count))
    duplicate_count <- as.integer(
        floor(feature_count * as.numeric(contract$duplicate_policy$fraction))
    )
    if (duplicate_count > 0L) {
        target <- seq.int(feature_count - duplicate_count + 1L, feature_count)
        feature_ids[target] <- feature_ids[seq_len(duplicate_count)]
    }
    internal_stride <- max(
        1L,
        as.integer(round(1 / as.numeric(parameters$internal_standard_fraction)))
    )
    internal <- seq_len(feature_count) %% internal_stride == 0L
    annotations <- sprintf("Synthetic metabolite %07d", seq_len(feature_count))
    unusual <- seq.int(13L, feature_count, by = 97L)
    annotations[unusual] <- paste0("Synthetic 'metabolite/", unusual, "'")
    data <- data.frame(
        assay = assignment,
        feature_id = feature_ids,
        annotation = annotations,
        internal_standard = internal,
        missingness_mechanism = mechanisms,
        mz = round(seq(75, 975, length.out = feature_count), 5L),
        retention_time = round(seq(0.1, 32, length.out = feature_count), 4L),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    sample_values <- as.data.frame(values, check.names = FALSE)
    names(sample_values) <- design$sample_id
    list(
        data = cbind(data, sample_values),
        design = design,
        effects = effects,
        missing_counts = missing$counts,
        duplicate_count = duplicate_count,
        assays = assays
    )
}

metabWorkloadWriteTable <- function(data, path) {
    utils::write.table(
        data,
        path,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = TRUE,
        na = "NA"
    )
    invisible(path)
}

metabWorkloadTruth <- function(contract, parameters, generated) {
    sample_columns <- generated$design$sample_id
    quantities <- as.matrix(generated$data[, sample_columns, drop = FALSE])
    list(
        schema = "multischolar.metabolomics_custom_truth",
        schema_version = "1.0.0",
        generator_version = parameters$generator_version,
        profile = parameters$profile,
        feature_rows = nrow(generated$data),
        sample_count = length(sample_columns),
        assays = as.list(generated$assays),
        sample_ids = as.list(sample_columns),
        group_assignments = as.list(generated$design$group),
        batch_assignments = as.list(generated$design$batch),
        internal_standard_count = sum(generated$data$internal_standard),
        duplicate_id_count = generated$duplicate_count,
        zero_count = sum(quantities == 0, na.rm = TRUE),
        missing_counts = as.list(generated$missing_counts),
        nonmissing_count = sum(!is.na(quantities)),
        finite_nonmissing = all(is.finite(quantities[!is.na(quantities)])),
        effect_count_per_direction = generated$effects$count,
        effect_log2 = as.numeric(parameters$effect_log2),
        column_names = as.list(names(generated$data)),
        privacy_classification = contract$privacy$classification,
        scale_mapping_digest = metabWorkloadDigest(list(
            dimensions = contract$dimensions,
            assay_mix = contract$assay_mix,
            scale_metadata = contract$privacy$scale_metadata
        ))
    )
}

metabWorkloadPrepare <- function(context) {
    contract <- context$contract
    parameters <- metabWorkloadParameters(contract)
    generated <- metabWorkloadPayload(contract, parameters)
    staging <- file.path(context$run_dir, "staging")
    dir.create(staging, recursive = TRUE, showWarnings = FALSE)
    payload_path <- file.path(staging, "synthetic-metabolomics.tsv")
    truth_path <- file.path(staging, "synthetic-metabolomics-truth.json")
    metabWorkloadWriteTable(generated$data, payload_path)
    jsonlite::write_json(
        metabWorkloadTruth(contract, parameters, generated),
        truth_path,
        auto_unbox = TRUE,
        pretty = TRUE,
        null = "null",
        na = "null",
        digits = 17
    )
    list(
        payload_path = payload_path,
        truth_path = truth_path,
        metadata = list(
            generated = TRUE,
            evidence_class = "generated_scaling",
            profile = parameters$profile,
            generator_version = parameters$generator_version,
            privacy_classification = contract$privacy$classification,
            feature_count = contract$dimensions$feature_count,
            sample_count = contract$dimensions$sample_count,
            assay_count = contract$dimensions$assay_count
        )
    )
}

metabWorkloadRead <- function(path) {
    utils::read.delim(
        path,
        check.names = FALSE,
        stringsAsFactors = FALSE,
        na.strings = "NA"
    )
}

metabWorkloadValidatePayload <- function(data, truth) {
    required <- c(
        "assay", "feature_id", "annotation", "internal_standard",
        "missingness_mechanism", "mz", "retention_time"
    )
    sample_columns <- unlist(truth$sample_ids, use.names = FALSE)
    valid <- is.data.frame(data) && all(required %in% names(data)) &&
        identical(nrow(data), as.integer(truth$feature_rows)) &&
        identical(names(data), unlist(truth$column_names, use.names = FALSE)) &&
        identical(sort(unique(data$assay)), sort(unlist(truth$assays))) &&
        all(sample_columns %in% names(data)) &&
        all(vapply(data[sample_columns], is.numeric, logical(1))) &&
        all(data$missingness_mechanism %in% c(
            "MCAR", "MAR", "MNAR", "left_censored"
        ))
    if (!valid) metabWorkloadAbort("generated metabolomics payload is invalid")
    quantities <- as.matrix(data[, sample_columns, drop = FALSE])
    finite <- all(is.finite(quantities[!is.na(quantities)]))
    duplicates <- sum(duplicated(data$feature_id))
    missing <- vapply(c("MCAR", "MAR", "MNAR", "left_censored"), \(name) {
        rows <- data$missingness_mechanism == name
        sum(is.na(quantities[rows, , drop = FALSE]))
    }, integer(1))
    identical(finite, isTRUE(truth$finite_nonmissing)) &&
        identical(duplicates, as.integer(truth$duplicate_id_count)) &&
        identical(unname(missing), as.integer(unlist(truth$missing_counts))) &&
        identical(sum(quantities == 0, na.rm = TRUE), as.integer(truth$zero_count))
}

metabWorkloadQueryEvidence <- function(data, truth, parameters) {
    count <- min(
        as.integer(parameters$bounded_query_features),
        nrow(data)
    )
    ids <- unique(data$feature_id)[seq_len(count)]
    timings <- numeric(as.integer(parameters$query_repetitions))
    rows <- 0L
    for (index in seq_along(timings)) {
        started <- proc.time()[["elapsed"]]
        selected <- data[data$feature_id %in% ids, , drop = FALSE]
        rows <- nrow(selected)
        timings[[index]] <- proc.time()[["elapsed"]] - started +
            .Machine$double.eps
    }
    list(list(
        query_id = "bounded_feature_lookup",
        rows = rows,
        p95_seconds = as.numeric(stats::quantile(
            timings,
            0.95,
            names = FALSE,
            type = 8
        ))
    ))
}

metabWorkloadRun <- function(context) {
    contract <- context$contract
    parameters <- metabWorkloadParameters(contract)
    truth <- jsonlite::read_json(context$truth_path, simplifyVector = TRUE)
    data <- metabWorkloadRead(context$payload_path)
    truth_valid <- metabWorkloadValidatePayload(data, truth)
    if (!isTRUE(truth_valid)) {
        metabWorkloadAbort("metabolomics workload differs from declared truth")
    }
    list(
        status = "passed",
        workflow_evidence = list(
            truth_valid = TRUE,
            evidence_class = "generated_scaling",
            profile = parameters$profile,
            generator_version = parameters$generator_version,
            assays = truth$assays,
            feature_rows = truth$feature_rows,
            sample_count = truth$sample_count,
            missing_counts = truth$missing_counts,
            duplicate_id_count = truth$duplicate_id_count,
            internal_standard_count = truth$internal_standard_count,
            scale_mapping_digest = truth$scale_mapping_digest,
            claim_boundary = paste(
                "custom schema workload/resource evidence only;",
                "not vendor or biological validation"
            )
        ),
        query_evidence = metabWorkloadQueryEvidence(data, truth, parameters),
        session_evidence = list(status = "memory_only", resource_count = 0L),
        report_evidence = list(
            status = "covered_by_separate_e2e_lane",
            file_count = 0L
        ),
        retained = data
    )
}

omicsWorkloadAdapter <- function() {
    list(
        adapter_id = .METAB_WORKLOAD_ADAPTER_ID,
        adapter_version = .METAB_WORKLOAD_ADAPTER_VERSION,
        supported_omics = "metabolomics",
        prepare = metabWorkloadPrepare,
        run = metabWorkloadRun
    )
}
