.LIPID_WORKLOAD_ADAPTER_ID <- "multischolar.lipidomics.v1"
.LIPID_WORKLOAD_ADAPTER_VERSION <- "1.1.0"

lipidWorkloadAbort <- function(message) stop(message, call. = FALSE)

lipidWorkloadDigest <- function(value) {
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

lipidWorkloadParameters <- function(contract) {
    parameters <- contract$adapter_parameters
    generated_required <- c(
        "source_kind",
        "generator_version", "profile", "internal_standard_fraction",
        "effect_fraction", "effect_log2", "batch_log2", "heteroscedasticity",
        "bounded_query_features", "query_repetitions", "scale_mapping_policy",
        "source_report_column_count", "source_byte_class"
    )
    fixture_required <- c(
        "source_kind", "profile", "bounded_query_features",
        "query_repetitions", "fixture_relpath", "oracle"
    )
    if (!is.list(parameters) || !is.character(parameters$source_kind) ||
        length(parameters$source_kind) != 1L) {
        lipidWorkloadAbort("lipidomics adapter parameters are unsupported")
    }
    if (identical(parameters$source_kind, "generated_custom")) {
        if (!setequal(names(parameters), generated_required) ||
            !identical(parameters$generator_version, "1.0.0")) {
            lipidWorkloadAbort("generated lipidomics parameters are unsupported")
        }
    } else if (identical(parameters$source_kind, "reviewed_lipidsearch_fixture")) {
        if (!setequal(names(parameters), fixture_required)) {
            lipidWorkloadAbort("fixture lipidomics parameters are unsupported")
        }
        lipidWorkloadValidateFixtureContract(contract, parameters)
    } else {
        lipidWorkloadAbort("lipidomics source kind is unsupported")
    }
    parameters
}

lipidWorkloadFixtureCatalogue <- function() {
    list(
        "tests/testdata/e2e/lipid_canonical/lipidsearch_lcms_pos.txt" =
            "LCMS_Pos",
        "tests/testdata/e2e/lipid_canonical/lipidsearch_lcms_neg.txt" =
            "LCMS_Neg",
        "tests/testdata/e2e/lipid_canonical/lipidsearch_gcms.txt" =
            "GCMS"
    )
}

lipidWorkloadValidateFixtureContract <- function(contract, parameters) {
    catalogue <- lipidWorkloadFixtureCatalogue()
    path <- parameters$fixture_relpath
    oracle_required <- c(
        "assay", "rows", "sample_columns", "quantity_sum", "first_lipid",
        "first_class", "first_adduct"
    )
    valid <- is.character(path) && length(path) == 1L &&
        path %in% names(catalogue) &&
        identical(contract$capability$input_format, "lipidsearch") &&
        identical(names(contract$assay_mix), unname(catalogue[[path]])) &&
        is.list(parameters$oracle) &&
        setequal(names(parameters$oracle), oracle_required) &&
        identical(parameters$oracle$assay, unname(catalogue[[path]])) &&
        identical(
            as.integer(contract$dimensions$feature_count),
            as.integer(parameters$oracle$rows)
        ) &&
        identical(
            as.integer(contract$dimensions$sample_count),
            length(unlist(parameters$oracle$sample_columns, use.names = FALSE))
        )
    if (!isTRUE(valid)) {
        lipidWorkloadAbort("reviewed LipidSearch fixture contract is inconsistent")
    }
    invisible(parameters)
}

lipidWorkloadResolveFixture <- function(relative_path) {
    roots <- normalizePath(
        c(
            getwd(),
            file.path(getwd(), ".."),
            file.path(getwd(), "..", ".."),
            file.path(getwd(), "..", "..", "..")
        ),
        mustWork = FALSE
    )
    candidates <- file.path(roots, relative_path)
    matches <- which(file.exists(candidates) & !dir.exists(candidates))
    if (!length(matches)) {
        lipidWorkloadAbort("allowlisted LipidSearch fixture could not be resolved")
    }
    root <- roots[[matches[[1L]]]]
    list(
        source = normalizePath(candidates[[matches[[1L]]]], mustWork = TRUE),
        fixture_root = normalizePath(
            file.path(root, "tests", "testdata", "e2e", "lipid_canonical"),
            mustWork = TRUE
        )
    )
}

lipidWorkloadAssays <- function(contract) {
    assays <- names(contract$assay_mix)
    supported <- c("LCMS_Pos", "LCMS_Neg", "GCMS")
    if (!length(assays) || !all(assays %in% supported) ||
        length(assays) != as.integer(contract$dimensions$assay_count)) {
        lipidWorkloadAbort("lipidomics assay mix is unsupported")
    }
    assays
}

lipidWorkloadDesign <- function(sample_count) {
    if (sample_count < 6L || sample_count %% 2L != 0L) {
        lipidWorkloadAbort("lipidomics sample_count must be even and at least six")
    }
    data.frame(
        sample_id = sprintf("LIPID_S%02d", seq_len(sample_count)),
        group = rep(c("CTRL", "TREAT"), each = sample_count / 2L),
        batch = rep(c("B1", "B2"), length.out = sample_count),
        stringsAsFactors = FALSE
    )
}

lipidWorkloadChemistry <- function(assignment, feature_index) {
    definitions <- list(
        LCMS_Pos = list(
            classes = c("PC", "PE", "TG"),
            adducts = c("[M+H]+", "[M+H]+", "[M+NH4]+"),
            ion_mode = "positive"
        ),
        LCMS_Neg = list(
            classes = c("PI", "PS", "FA"),
            adducts = c("[M-H]-", "[M-H]-", "[M-H]-"),
            ion_mode = "negative"
        ),
        GCMS = list(
            classes = c("FA", "FAME", "Sterol"),
            adducts = c("[M]+.", "[M]+.", "[M]+."),
            ion_mode = "EI"
        )
    )
    rows <- lapply(seq_along(assignment), \(index) {
        definition <- definitions[[assignment[[index]]]]
        chemistry_index <- (feature_index[[index]] - 1L) %%
            length(definition$classes) + 1L
        list(
            lipid_class = definition$classes[[chemistry_index]],
            adduct = definition$adducts[[chemistry_index]],
            ion_mode = definition$ion_mode
        )
    })
    data.frame(
        lipid_class = vapply(rows, `[[`, character(1), "lipid_class"),
        adduct = vapply(rows, `[[`, character(1), "adduct"),
        ion_mode = vapply(rows, `[[`, character(1), "ion_mode"),
        stringsAsFactors = FALSE
    )
}

lipidWorkloadIntensity <- function(contract, parameters, assignment, design) {
    means <- unlist(contract$distributions$assay_log2_means, use.names = TRUE)
    feature_count <- length(assignment)
    baseline <- stats::rnorm(
        feature_count,
        means[assignment],
        as.numeric(contract$distributions$base_log2_sd)
    )
    effect_count <- max(1L, as.integer(floor(
        feature_count * as.numeric(parameters$effect_fraction)
    )))
    effects <- numeric(feature_count)
    effects[seq_len(effect_count)] <- as.numeric(parameters$effect_log2)
    effects[seq.int(effect_count + 1L, 2L * effect_count)] <-
        -as.numeric(parameters$effect_log2)
    treated <- as.integer(design$group == "TREAT")
    batch <- ifelse(design$batch == "B2", as.numeric(parameters$batch_log2), 0)
    hetero <- 1 + as.numeric(parameters$heteroscedasticity) *
        abs(baseline - mean(baseline)) / max(stats::sd(baseline), 1e-8)
    noise <- matrix(
        stats::rnorm(feature_count * nrow(design)),
        nrow = feature_count
    ) * as.numeric(contract$distributions$noise_log2_sd) * hetero
    values <- baseline + effects %o% treated
    values <- sweep(values, 2L, batch, "+") + noise
    list(
        values = round(pmax(2^values, 0), 6L),
        effects = effects,
        effect_count = effect_count
    )
}

lipidWorkloadMissingness <- function(values, mechanisms, assignment, design, contract) {
    counts <- stats::setNames(integer(4L), c(
        "MCAR", "MAR", "MNAR", "left_censored"
    ))
    thresholds <- unlist(
        contract$censoring$assay_intensity_thresholds,
        use.names = TRUE
    )
    for (row in seq_len(nrow(values))) {
        mechanism <- mechanisms[[row]]
        columns <- switch(mechanism,
            MCAR = which((row + seq_len(ncol(values))) %% 19L == 0L),
            MAR = which(design$batch == "B2" &
                (row + seq_len(ncol(values))) %% 7L == 0L),
            MNAR = order(values[row, ], method = "radix")[[1L]],
            left_censored = which(values[row, ] < thresholds[[assignment[[row]]]])
        )
        if (length(columns)) {
            values[row, columns] <- NA_real_
            counts[[mechanism]] <- counts[[mechanism]] + length(columns)
        }
    }
    list(values = values, counts = counts)
}

lipidWorkloadPayload <- function(contract, parameters) {
    feature_count <- as.integer(contract$dimensions$feature_count)
    assays <- lipidWorkloadAssays(contract)
    design <- lipidWorkloadDesign(as.integer(contract$dimensions$sample_count))
    assignment <- rep(assays, length.out = feature_count)
    feature_index <- stats::ave(
        seq_len(feature_count),
        assignment,
        FUN = seq_along
    )
    chemistry <- lipidWorkloadChemistry(assignment, feature_index)
    intensity <- lipidWorkloadIntensity(
        contract,
        parameters,
        assignment,
        design
    )
    mechanisms <- rep(c("MCAR", "MAR", "MNAR", "left_censored"),
        length.out = feature_count
    )
    missing <- lipidWorkloadMissingness(
        intensity$values,
        mechanisms,
        assignment,
        design,
        contract
    )
    values <- missing$values
    values[[1L, 1L]] <- 0
    lipid_ids <- sprintf(
        "%s %02d:%d",
        chemistry$lipid_class,
        30L + feature_index %% 25L,
        feature_index %% 7L
    )
    duplicate_count <- as.integer(floor(
        feature_count * as.numeric(contract$duplicate_policy$fraction)
    ))
    if (duplicate_count > 0L) {
        target <- seq.int(feature_count - duplicate_count + 1L, feature_count)
        lipid_ids[target] <- lipid_ids[seq_len(duplicate_count)]
    }
    stride <- max(1L, as.integer(round(
        1 / as.numeric(parameters$internal_standard_fraction)
    )))
    internal <- seq_len(feature_count) %% stride == 0L
    lipid_ids[internal] <- paste0("IS_", lipid_ids[internal], "_d7")
    data <- data.frame(
        assay = assignment,
        lipid_id = lipid_ids,
        lipid_class = chemistry$lipid_class,
        adduct = chemistry$adduct,
        ion_mode = chemistry$ion_mode,
        internal_standard = internal,
        missingness_mechanism = mechanisms,
        mz = round(seq(150, 1100, length.out = feature_count), 5L),
        retention_time = round(seq(0.1, 35, length.out = feature_count), 4L),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    sample_values <- as.data.frame(values, check.names = FALSE)
    names(sample_values) <- design$sample_id
    list(
        data = cbind(data, sample_values),
        design = design,
        assays = assays,
        missing_counts = missing$counts,
        duplicate_count = duplicate_count,
        effect_count = intensity$effect_count
    )
}

lipidWorkloadTruth <- function(contract, parameters, generated) {
    samples <- generated$design$sample_id
    quantities <- as.matrix(generated$data[samples])
    list(
        schema = "multischolar.lipidomics_custom_truth",
        schema_version = "1.0.0",
        generator_version = parameters$generator_version,
        profile = parameters$profile,
        feature_rows = nrow(generated$data),
        sample_count = length(samples),
        assays = as.list(generated$assays),
        sample_ids = as.list(samples),
        group_assignments = as.list(generated$design$group),
        batch_assignments = as.list(generated$design$batch),
        lipid_classes = as.list(sort(
            unique(generated$data$lipid_class),
            method = "radix"
        )),
        adducts = as.list(sort(
            unique(generated$data$adduct),
            method = "radix"
        )),
        ion_modes = as.list(sort(
            unique(generated$data$ion_mode),
            method = "radix"
        )),
        internal_standard_count = sum(generated$data$internal_standard),
        duplicate_id_count = generated$duplicate_count,
        zero_count = sum(quantities == 0, na.rm = TRUE),
        missing_counts = as.list(generated$missing_counts),
        finite_nonmissing = all(is.finite(quantities[!is.na(quantities)])),
        effect_count_per_direction = generated$effect_count,
        effect_log2 = as.numeric(parameters$effect_log2),
        column_names = as.list(names(generated$data)),
        privacy_classification = contract$privacy$classification,
        scale_mapping_digest = digest::digest(
            paste(
                contract$dimensions$feature_count,
                contract$dimensions$sample_count,
                paste(names(contract$assay_mix), collapse = ","),
                contract$privacy$classification,
                parameters$scale_mapping_policy,
                parameters$source_report_column_count,
                parameters$source_byte_class,
                jsonlite::toJSON(
                    contract$privacy$scale_metadata,
                    auto_unbox = TRUE,
                    null = "null",
                    digits = 17
                ),
                sep = "|"
            ),
            algo = "sha256",
            serialize = FALSE
        )
    )
}

lipidWorkloadPrepareGenerated <- function(context, parameters) {
    contract <- context$contract
    generated <- lipidWorkloadPayload(contract, parameters)
    staging <- file.path(context$run_dir, "staging")
    dir.create(staging, recursive = TRUE, showWarnings = FALSE)
    payload_path <- file.path(staging, "synthetic-lipidomics.tsv")
    truth_path <- file.path(staging, "synthetic-lipidomics-truth.json")
    utils::write.table(
        generated$data,
        payload_path,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        na = "NA"
    )
    jsonlite::write_json(
        lipidWorkloadTruth(contract, parameters, generated),
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

lipidWorkloadPrepareFixture <- function(context, parameters) {
    staging <- file.path(context$run_dir, "staging")
    dir.create(staging, recursive = TRUE, showWarnings = FALSE)
    resolved <- lipidWorkloadResolveFixture(parameters$fixture_relpath)
    source <- resolved$source
    fixture_root <- resolved$fixture_root
    if (!startsWith(source, paste0(fixture_root, .Platform$file.sep))) {
        lipidWorkloadAbort("LipidSearch fixture resolved outside its allowlisted root")
    }
    payload_path <- file.path(staging, basename(source))
    if (!file.copy(source, payload_path, overwrite = FALSE)) {
        lipidWorkloadAbort("failed to stage the reviewed LipidSearch fixture")
    }
    truth_path <- file.path(staging, "lipidsearch-fixture-truth.json")
    truth <- c(
        list(
            schema = "multischolar.lipidomics_fixture_truth",
            schema_version = "1.0.0",
            evidence_class = "independent_fixture",
            profile = parameters$profile
        ),
        parameters$oracle
    )
    jsonlite::write_json(
        truth,
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
            generated = FALSE,
            evidence_class = "independent_fixture",
            profile = parameters$profile,
            input_format = "lipidsearch",
            assay = parameters$oracle$assay,
            feature_count = context$contract$dimensions$feature_count,
            sample_count = context$contract$dimensions$sample_count
        )
    )
}

lipidWorkloadPrepare <- function(context) {
    parameters <- lipidWorkloadParameters(context$contract)
    if (identical(parameters$source_kind, "reviewed_lipidsearch_fixture")) {
        return(lipidWorkloadPrepareFixture(context, parameters))
    }
    lipidWorkloadPrepareGenerated(context, parameters)
}

lipidWorkloadQueryEvidence <- function(data, id_column, parameters, query_id) {
    count <- min(
        as.integer(parameters$bounded_query_features),
        length(unique(data[[id_column]]))
    )
    ids <- unique(data[[id_column]])[seq_len(count)]
    timings <- numeric(as.integer(parameters$query_repetitions))
    rows <- 0L
    for (index in seq_along(timings)) {
        started <- proc.time()[["elapsed"]]
        selected <- data[data[[id_column]] %in% ids, , drop = FALSE]
        rows <- nrow(selected)
        timings[[index]] <- proc.time()[["elapsed"]] - started +
            .Machine$double.eps
    }
    list(list(
        query_id = query_id,
        rows = rows,
        p95_seconds = as.numeric(stats::quantile(
            timings,
            0.95,
            names = FALSE,
            type = 8
        ))
    ))
}

lipidWorkloadRunGenerated <- function(context, parameters) {
    contract <- context$contract
    truth <- jsonlite::read_json(context$truth_path, simplifyVector = TRUE)
    data <- utils::read.delim(
        context$payload_path,
        check.names = FALSE,
        stringsAsFactors = FALSE,
        na.strings = "NA"
    )
    required <- c(
        "assay", "lipid_id", "lipid_class", "adduct", "ion_mode",
        "internal_standard", "missingness_mechanism", "mz", "retention_time"
    )
    valid <- all(required %in% names(data)) &&
        identical(nrow(data), as.integer(truth$feature_rows)) &&
        identical(names(data), unlist(truth$column_names, use.names = FALSE)) &&
        identical(sort(unique(data$lipid_class)), sort(truth$lipid_classes)) &&
        identical(sort(unique(data$adduct)), sort(truth$adducts)) &&
        identical(sort(unique(data$ion_mode)), sort(truth$ion_modes))
    if (!isTRUE(valid)) lipidWorkloadAbort("generated lipidomics truth mismatch")
    list(
        status = "passed",
        workflow_evidence = list(
            truth_valid = TRUE,
            evidence_class = "generated_scaling",
            profile = parameters$profile,
            assays = truth$assays,
            lipid_classes = truth$lipid_classes,
            adducts = truth$adducts,
            ion_modes = truth$ion_modes,
            missing_counts = truth$missing_counts,
            scale_mapping_digest = truth$scale_mapping_digest,
            claim_boundary = paste(
                "custom schema workload/resource evidence only;",
                "not vendor or biological validation"
            )
        ),
        query_evidence = lipidWorkloadQueryEvidence(
            data,
            "lipid_id",
            parameters,
            "bounded_lipid_lookup"
        ),
        session_evidence = list(status = "memory_only", resource_count = 0L),
        report_evidence = list(status = "separate_e2e", file_count = 0L),
        retained = data
    )
}

lipidWorkloadRunFixture <- function(context, parameters) {
    if (is.null(context$package_library)) {
        loadNamespace("MultiScholaR")
    } else {
        loadNamespace("MultiScholaR", lib.loc = context$package_library)
    }
    importer <- getExportedValue("MultiScholaR", "importLipidSearchData")
    imported <- suppressMessages(importer(context$payload_path))
    truth <- jsonlite::read_json(context$truth_path, simplifyVector = TRUE)
    sample_columns <- unlist(truth$sample_columns, use.names = FALSE)
    data <- imported$data
    values <- as.matrix(data[sample_columns])
    valid <- identical(imported$format, "lipidsearch") &&
        identical(nrow(data), as.integer(truth$rows)) &&
        identical(imported$sample_columns, sample_columns) &&
        isTRUE(all.equal(sum(values), as.numeric(truth$quantity_sum))) &&
        identical(data$LipidName[[1L]], truth$first_lipid) &&
        identical(data$LipidClass[[1L]], truth$first_class) &&
        identical(data$IonType[[1L]], truth$first_adduct)
    if (!isTRUE(valid)) {
        lipidWorkloadAbort("reviewed LipidSearch fixture truth mismatch")
    }
    list(
        status = "passed",
        workflow_evidence = list(
            truth_valid = TRUE,
            evidence_class = "independent_fixture",
            profile = parameters$profile,
            input_format = "lipidsearch",
            assay = truth$assay,
            rows = truth$rows,
            sample_count = length(sample_columns),
            quantity_sum = truth$quantity_sum,
            claim_boundary = paste(
                "direct LipidSearch reader and resource evidence only;",
                "not representative-scale or biological validation"
            )
        ),
        query_evidence = lipidWorkloadQueryEvidence(
            data,
            "LipidName",
            parameters,
            "bounded_lipidsearch_lookup"
        ),
        session_evidence = list(status = "memory_only", resource_count = 0L),
        report_evidence = list(status = "separate_e2e", file_count = 0L),
        retained = data
    )
}

lipidWorkloadRun <- function(context) {
    parameters <- lipidWorkloadParameters(context$contract)
    if (identical(parameters$source_kind, "reviewed_lipidsearch_fixture")) {
        return(lipidWorkloadRunFixture(context, parameters))
    }
    lipidWorkloadRunGenerated(context, parameters)
}

omicsWorkloadAdapter <- function() {
    list(
        adapter_id = .LIPID_WORKLOAD_ADAPTER_ID,
        adapter_version = .LIPID_WORKLOAD_ADAPTER_VERSION,
        supported_omics = "lipidomics",
        prepare = lipidWorkloadPrepare,
        run = lipidWorkloadRun
    )
}
