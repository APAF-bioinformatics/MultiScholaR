.PROTEOMICS_WORKLOAD_ADAPTER_ID <- "multischolar.proteomics.nondia.v1"
.PROTEOMICS_WORKLOAD_ADAPTER_VERSION <- "1.0.0"

proteomicsWorkloadAbort <- function(message) {
    stop(message, call. = FALSE)
}

proteomicsWorkloadScalar <- function(value, label) {
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
        proteomicsWorkloadAbort(sprintf("%s must be one finite number", label))
    }
    as.numeric(value)
}

proteomicsWorkloadParameters <- function(contract) {
    parameters <- contract$adapter_parameters
    required <- c(
        "format",
        "effect_fraction",
        "effect_log2",
        "base_log2_mean",
        "base_log2_sd",
        "noise_log2_sd",
        "bounded_query_features",
        "query_repetitions",
        "e2e_lane"
    )
    if (!is.list(parameters) || !setequal(names(parameters), required)) {
        proteomicsWorkloadAbort(
            "proteomics adapter parameters do not match the versioned contract"
        )
    }
    supported <- c("maxquant", "fragpipe", "pd_tmt")
    if (!parameters$format %in% supported) {
        proteomicsWorkloadAbort(sprintf(
            "unsupported proteomics workload format '%s'",
            as.character(parameters$format)
        ))
    }
    parameters
}

proteomicsWorkloadChannels <- function(sample_count) {
    channels <- c(
        "126", "127N", "127C", "128N", "128C", "129N", "129C",
        "130N", "130C", "131N", "131C", "132N", "132C", "133N",
        "133C", "134N", "134C", "135N"
    )
    if (sample_count > length(channels)) {
        proteomicsWorkloadAbort("pd_tmt workloads support at most 18 samples")
    }
    channels[seq_len(sample_count)]
}

proteomicsWorkloadSamples <- function(sample_count) {
    if (sample_count < 4L || sample_count %% 2L != 0L) {
        proteomicsWorkloadAbort("sample_count must be an even integer of at least four")
    }
    group_size <- sample_count / 2L
    c(
        sprintf("CTRL_%02d", seq_len(group_size)),
        sprintf("TREAT_%02d", seq_len(group_size))
    )
}

proteomicsWorkloadQuantities <- function(contract, parameters) {
    feature_count <- as.integer(contract$dimensions$feature_count)
    sample_count <- as.integer(contract$dimensions$sample_count)
    effect_fraction <- proteomicsWorkloadScalar(
        parameters$effect_fraction,
        "effect_fraction"
    )
    if (effect_fraction <= 0 || effect_fraction >= 0.5) {
        proteomicsWorkloadAbort("effect_fraction must be between zero and 0.5")
    }

    effect_count <- max(1L, as.integer(floor(feature_count * effect_fraction)))
    effects <- numeric(feature_count)
    effects[seq_len(effect_count)] <- parameters$effect_log2
    down_index <- seq.int(effect_count + 1L, 2L * effect_count)
    effects[down_index] <- -parameters$effect_log2

    baseline <- stats::rnorm(
        feature_count,
        mean = parameters$base_log2_mean,
        sd = parameters$base_log2_sd
    )
    sample_offsets <- stats::rnorm(sample_count, mean = 0, sd = 0.08)
    noise <- matrix(
        stats::rnorm(
            feature_count * sample_count,
            mean = 0,
            sd = parameters$noise_log2_sd
        ),
        nrow = feature_count,
        ncol = sample_count
    )
    treated <- c(rep(0, sample_count / 2L), rep(1, sample_count / 2L))
    log2_values <- baseline + effects %o% treated
    log2_values <- sweep(log2_values, 2L, sample_offsets, "+") + noise
    quantities <- round(pmax(2^log2_values, 0), digits = 6L)

    list(
        values = quantities,
        effects = effects,
        effect_count = effect_count
    )
}

proteomicsWorkloadWriteTable <- function(data, path) {
    numeric_columns <- vapply(data, is.numeric, logical(1))
    data[numeric_columns] <- lapply(data[numeric_columns], \(column) {
        formatC(column, format = "f", digits = 6L)
    })
    utils::write.table(
        data,
        file = path,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = TRUE,
        na = ""
    )
    invisible(path)
}

proteomicsWorkloadMaxQuant <- function(feature_ids, gene_ids, samples, values) {
    data <- data.frame(
        Protein.IDs = feature_ids,
        Gene.names = gene_ids,
        Potential.contaminant = "",
        Reverse = "",
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    quantities <- as.data.frame(values, check.names = FALSE)
    names(quantities) <- paste0("LFQ.intensity.", samples)
    cbind(data, quantities)
}

proteomicsWorkloadFragPipe <- function(feature_ids, gene_ids, samples, values) {
    data <- data.frame(
        `Protein ID` = feature_ids,
        Gene = gene_ids,
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    quantities <- as.data.frame(values, check.names = FALSE)
    names(quantities) <- paste(samples, "MaxLFQ Intensity")
    cbind(data, quantities)
}

proteomicsWorkloadPdTmt <- function(feature_ids, samples, values) {
    data <- data.frame(
        Annotated.Sequence = sprintf("PEPTIDE%06dK", seq_along(feature_ids)),
        Master.Protein.Accessions = feature_ids,
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    channels <- proteomicsWorkloadChannels(length(samples))
    quantities <- as.data.frame(values, check.names = FALSE)
    names(quantities) <- sprintf(
        "Abundance: F%d: %s, %s",
        seq_along(samples),
        channels,
        samples
    )
    cbind(data, quantities)
}

proteomicsWorkloadPayload <- function(format, feature_ids, gene_ids, samples, values) {
    switch(
        format,
        maxquant = proteomicsWorkloadMaxQuant(
            feature_ids,
            gene_ids,
            samples,
            values
        ),
        fragpipe = proteomicsWorkloadFragPipe(
            feature_ids,
            gene_ids,
            samples,
            values
        ),
        pd_tmt = proteomicsWorkloadPdTmt(feature_ids, samples, values)
    )
}

proteomicsWorkloadPrepare <- function(context) {
    contract <- context$contract
    parameters <- proteomicsWorkloadParameters(contract)
    feature_count <- as.integer(contract$dimensions$feature_count)
    sample_count <- as.integer(contract$dimensions$sample_count)
    samples <- proteomicsWorkloadSamples(sample_count)
    feature_ids <- sprintf("P%06d", seq_len(feature_count))
    gene_ids <- sprintf("GENE%06d", seq_len(feature_count))
    generated <- proteomicsWorkloadQuantities(contract, parameters)
    payload <- proteomicsWorkloadPayload(
        parameters$format,
        feature_ids,
        gene_ids,
        samples,
        generated$values
    )

    staging <- file.path(context$run_dir, "staging")
    dir.create(staging, recursive = TRUE, showWarnings = FALSE)
    extension <- if (identical(parameters$format, "maxquant")) ".txt" else ".tsv"
    payload_path <- file.path(staging, paste0("synthetic-proteomics", extension))
    truth_path <- file.path(staging, "synthetic-proteomics-truth.json")
    proteomicsWorkloadWriteTable(payload, payload_path)
    truth <- list(
        schema = "multischolar.proteomics_nondia_truth",
        schema_version = "1.0.0",
        format = parameters$format,
        feature_count = feature_count,
        sample_count = sample_count,
        long_row_count = feature_count * sample_count,
        quantity_sum = sum(generated$values),
        quantity_min = min(generated$values),
        quantity_max = max(generated$values),
        effect_count_per_direction = generated$effect_count,
        effect_log2 = parameters$effect_log2,
        sample_groups = list(control = sample_count / 2L, treatment = sample_count / 2L)
    )
    jsonlite::write_json(
        truth,
        truth_path,
        auto_unbox = TRUE,
        pretty = TRUE,
        digits = 17
    )
    list(
        payload_path = payload_path,
        truth_path = truth_path,
        metadata = list(
            generated = TRUE,
            evidence_class = "generated_scaling",
            format = parameters$format,
            feature_count = feature_count,
            sample_count = sample_count
        )
    )
}

proteomicsWorkloadImporter <- function(format) {
    namespace <- asNamespace("MultiScholaR")
    function_name <- switch(
        format,
        maxquant = "importMaxQuantData",
        fragpipe = "importFragPipeData",
        pd_tmt = "importProteomeDiscovererTMTData"
    )
    get(function_name, envir = namespace, inherits = FALSE)
}

proteomicsWorkloadImportSummary <- function(imported) {
    data <- as.data.frame(imported$data)
    mapping <- imported$column_mapping
    feature_ids <- as.character(data[[mapping$protein_col]])
    samples <- as.character(data[[mapping$run_col]])
    quantities <- as.numeric(data[[mapping$quantity_col]])
    list(
        data = data,
        feature_ids = feature_ids,
        samples = samples,
        quantities = quantities,
        evidence = list(
            data_type = imported$data_type,
            rows = nrow(data),
            columns = ncol(data),
            feature_count = length(unique(feature_ids)),
            sample_count = length(unique(samples)),
            column_names = as.list(names(data)),
            quantity_sum = sum(quantities),
            quantity_min = min(quantities),
            quantity_max = max(quantities),
            quantity_na_count = sum(is.na(quantities))
        )
    )
}

proteomicsWorkloadDirectionEvidence <- function(summary, truth) {
    log2_values <- log2(summary$quantities + 1)
    control <- grepl("CTRL", summary$samples, fixed = TRUE)
    treatment <- grepl("TREAT", summary$samples, fixed = TRUE)
    feature_order <- unique(summary$feature_ids)
    groupMean <- function(selected) {
        feature_ids <- summary$feature_ids[selected]
        sums <- rowsum(log2_values[selected], feature_ids, reorder = FALSE)
        counts <- rowsum(rep.int(1, length(feature_ids)), feature_ids, reorder = FALSE)
        means <- as.numeric(sums / counts)
        stats::setNames(means, rownames(sums))
    }
    control_means <- groupMean(control)
    treatment_means <- groupMean(treatment)
    effects <- treatment_means[feature_order] - control_means[feature_order]
    effect_count <- as.integer(truth$effect_count_per_direction)
    up_index <- seq_len(effect_count)
    down_index <- seq.int(effect_count + 1L, 2L * effect_count)
    unaffected_index <- seq.int(2L * effect_count + 1L, length(effects))
    margin <- max(0.25, abs(as.numeric(truth$effect_log2)) / 2)
    expected_directions <- all(effects[up_index] > margin) &&
        all(effects[down_index] < -margin)
    unaffected_bounded <- all(abs(effects[unaffected_index]) < margin)
    list(
        valid = expected_directions && unaffected_bounded,
        expected_up_count = effect_count,
        expected_down_count = effect_count,
        observed_up_count = sum(effects > margin),
        observed_down_count = sum(effects < -margin),
        unaffected_max_abs_log2fc = max(abs(effects[-c(up_index, down_index)]))
    )
}

proteomicsWorkloadTruthValid <- function(summary, truth, direction) {
    tolerance <- max(1, abs(as.numeric(truth$quantity_sum))) * 1e-12
    identical(summary$evidence$data_type, "protein") &&
        identical(summary$evidence$rows, as.integer(truth$long_row_count)) &&
        identical(summary$evidence$feature_count, as.integer(truth$feature_count)) &&
        identical(summary$evidence$sample_count, as.integer(truth$sample_count)) &&
        summary$evidence$quantity_na_count == 0L &&
        abs(summary$evidence$quantity_sum - truth$quantity_sum) <= tolerance &&
        isTRUE(direction$valid)
}

proteomicsWorkloadQueryEvidence <- function(summary, parameters) {
    query_features <- min(
        as.integer(parameters$bounded_query_features),
        summary$evidence$feature_count
    )
    repetitions <- as.integer(parameters$query_repetitions)
    selected_ids <- unique(summary$feature_ids)[seq_len(query_features)]
    timings <- numeric(repetitions)
    row_count <- 0L
    for (index in seq_len(repetitions)) {
        started <- proc.time()[["elapsed"]]
        selected <- summary$data[summary$feature_ids %in% selected_ids, , drop = FALSE]
        row_count <- nrow(selected)
        timings[[index]] <- proc.time()[["elapsed"]] - started + .Machine$double.eps
    }
    list(list(
        query_id = "bounded_feature_lookup",
        rows = row_count,
        p95_seconds = as.numeric(stats::quantile(
            timings,
            probs = 0.95,
            names = FALSE,
            type = 8
        ))
    ))
}

proteomicsWorkloadRun <- function(context) {
    parameters <- proteomicsWorkloadParameters(context$contract)
    truth <- jsonlite::read_json(context$truth_path, simplifyVector = TRUE)
    importer <- proteomicsWorkloadImporter(parameters$format)
    imported <- suppressWarnings(suppressMessages(importer(context$payload_path)))
    summary <- proteomicsWorkloadImportSummary(imported)
    direction <- proteomicsWorkloadDirectionEvidence(summary, truth)
    truth_valid <- proteomicsWorkloadTruthValid(summary, truth, direction)
    if (!identical(as.character(truth$format), parameters$format) || !truth_valid) {
        proteomicsWorkloadAbort(
            "imported proteomics workload does not match declared truth"
        )
    }

    list(
        status = "passed",
        workflow_evidence = list(
            truth_valid = truth_valid,
            evidence_class = "generated_scaling",
            format = parameters$format,
            import = summary$evidence,
            differential_direction = direction,
            reachable_stage_claim = "import_and_bounded_memory_query_only",
            e2e_lane = parameters$e2e_lane
        ),
        query_evidence = proteomicsWorkloadQueryEvidence(summary, parameters),
        session_evidence = list(status = "memory_only", resource_count = 0L),
        report_evidence = list(
            status = "covered_by_separate_browser_lane",
            file_count = 0L
        ),
        retained = imported$data
    )
}

omicsWorkloadAdapter <- function() {
    list(
        adapter_id = .PROTEOMICS_WORKLOAD_ADAPTER_ID,
        adapter_version = .PROTEOMICS_WORKLOAD_ADAPTER_VERSION,
        supported_omics = "proteomics",
        prepare = proteomicsWorkloadPrepare,
        run = proteomicsWorkloadRun
    )
}
