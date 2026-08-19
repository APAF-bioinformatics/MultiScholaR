.OMICS_PARITY_MANIFEST_SCHEMA <- "1.0.0"
.OMICS_PARITY_ORACLE_SCHEMA <- "1.0.0"
.OMICS_PARITY_RESULT_SCHEMA <- "1.0.0"

omicsParityRepoPath <- function(...) {
    file.path(normalizePath(testthat::test_path("..", ".."), mustWork = TRUE), ...)
}

omicsParityManifestPath <- function() {
    testthat::test_path("..", "testdata", "omics-parity", "scenarios.json")
}

omicsParityOraclePath <- function() {
    testthat::test_path("..", "testdata", "omics-parity", "dia-memory-oracle.json")
}

omicsParityReadJson <- function(path, expected_schema, document_name) {
    if (!file.exists(path)) {
        stop(sprintf("%s not found: %s", document_name, path), call. = FALSE)
    }
    document <- jsonlite::read_json(path, simplifyVector = FALSE)
    if (!identical(document$schema_version, expected_schema)) {
        observed_schema <- document$schema_version
        if (is.null(observed_schema)) {
            observed_schema <- "<missing>"
        }
        stop(
            sprintf(
                "Unsupported %s schema version '%s'; expected '%s'",
                document_name,
                observed_schema,
                expected_schema
            ),
            call. = FALSE
        )
    }
    document
}

omicsParityReadManifest <- function(path = omicsParityManifestPath()) {
    manifest <- omicsParityReadJson(
        path,
        .OMICS_PARITY_MANIFEST_SCHEMA,
        "omics parity manifest"
    )
    required <- c(
        "corpus_version", "capability_id", "backend", "oracle_path", "execution",
        "thread_controls",
        "comparison_contract", "release_gates", "disk_categories", "scenarios",
        "bounded_queries", "scientific_targets", "cross_omic_regression_targets"
    )
    missing <- setdiff(required, names(manifest))
    if (length(missing)) {
        stop(
            sprintf("Omics parity manifest is missing: %s", paste(missing, collapse = ", ")),
            call. = FALSE
        )
    }
    if (!identical(manifest$backend, "memory")) {
        stop("The baseline manifest must use the memory backend", call. = FALSE)
    }

    scenario_ids <- vapply(manifest$scenarios, `[[`, character(1), "scenario_id")
    target_ids <- vapply(manifest$scientific_targets, `[[`, character(1), "target_id")
    if (anyDuplicated(scenario_ids) || anyDuplicated(target_ids)) {
        stop("Scenario and scientific target identifiers must be unique", call. = FALSE)
    }
    manifest
}

omicsParityReadOracle <- function(path = omicsParityOraclePath()) {
    oracle <- omicsParityReadJson(
        path,
        .OMICS_PARITY_ORACLE_SCHEMA,
        "DIA-NN memory oracle"
    )
    oracle_ids <- vapply(oracle$scenarios, `[[`, character(1), "scenario_id")
    if (anyDuplicated(oracle_ids)) {
        stop("DIA-NN memory oracle scenario identifiers must be unique", call. = FALSE)
    }
    oracle
}

omicsParitySha256File <- function(path) {
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

omicsParityScalarText <- function(value) {
    if (length(value) != 1L) {
        stop("Canonical values must be scalar", call. = FALSE)
    }
    if (is.na(value)) {
        return(if (is.nan(value)) "<NaN>" else "<NA>")
    }
    if (is.numeric(value)) {
        if (is.infinite(value)) {
            return(if (value > 0) "<Inf>" else "<-Inf>")
        }
        return(sprintf("%.17g", value))
    }
    if (is.logical(value)) {
        return(if (value) "TRUE" else "FALSE")
    }
    enc2utf8(as.character(value))
}

omicsParityTableSha256 <- function(data, stable_keys, ordered_columns) {
    missing_keys <- setdiff(stable_keys, names(data))
    missing_columns <- setdiff(ordered_columns, names(data))
    if (length(missing_keys) || length(missing_columns)) {
        stop(
            sprintf(
                "Cannot canonicalize table; missing keys/columns: %s",
                paste(unique(c(missing_keys, missing_columns)), collapse = ", ")
            ),
            call. = FALSE
        )
    }

    order_values <- lapply(data[stable_keys], \(column) {
        vapply(column, omicsParityScalarText, character(1))
    })
    row_order <- do.call(order, c(order_values, list(na.last = TRUE, method = "radix")))
    canonical <- data[row_order, ordered_columns, drop = FALSE]
    rows <- apply(canonical, 1L, \(row) {
        paste(vapply(row, omicsParityScalarText, character(1)), collapse = "\t")
    })
    payload <- paste(c(paste(ordered_columns, collapse = "\t"), rows), collapse = "\n")
    digest::digest(payload, algo = "sha256", serialize = FALSE)
}

omicsParitySummarizeDiann <- function(import_result, fixture_path, contract) {
    data <- as.data.frame(import_result$data, stringsAsFactors = FALSE)
    quantity_column <- import_result$column_mapping$quantity_col
    quantity <- data[[quantity_column]]
    stable_keys <- unlist(contract$stable_keys, use.names = FALSE)
    ordered_columns <- unlist(contract$ordered_columns, use.names = FALSE)

    list(
        input_sha256 = omicsParitySha256File(fixture_path),
        output_sha256 = omicsParityTableSha256(data, stable_keys, ordered_columns),
        rows = nrow(data),
        columns = ncol(data),
        column_names = as.list(names(data)),
        column_classes = as.list(vapply(data, \(column) class(column)[[1L]], character(1))),
        runs = as.list(sort(unique(as.character(data$Run)))),
        proteins = as.list(sort(unique(as.character(data$Protein.Group)))),
        peptides = as.list(sort(unique(as.character(data$Stripped.Sequence)))),
        quantity_column = quantity_column,
        quantity_sum = sum(quantity, na.rm = TRUE),
        quantity_min = min(quantity, na.rm = TRUE),
        quantity_max = max(quantity, na.rm = TRUE),
        na_quantity = sum(is.na(quantity)),
        nan_quantity = sum(is.nan(quantity)),
        positive_infinite_quantity = sum(is.infinite(quantity) & quantity > 0, na.rm = TRUE),
        negative_infinite_quantity = sum(is.infinite(quantity) & quantity < 0, na.rm = TRUE)
    )
}

omicsParityCompareSummary <- function(actual, expected, contract) {
    exact_fields <- unlist(contract$exact_fields, use.names = FALSE)
    numeric_fields <- names(contract$numeric_fields)
    errors <- character()

    for (field in exact_fields) {
        if (!identical(actual[[field]], expected[[field]])) {
            errors <- c(errors, sprintf("exact field '%s' differs", field))
        }
    }
    for (field in numeric_fields) {
        tolerance <- contract$numeric_fields[[field]]
        difference <- abs(as.numeric(actual[[field]]) - as.numeric(expected[[field]]))
        allowed <- tolerance$absolute + tolerance$relative * abs(as.numeric(expected[[field]]))
        if (is.na(difference) || difference > allowed) {
            errors <- c(
                errors,
                sprintf(
                    "numeric field '%s' differs by %.17g (allowed %.17g)",
                    field,
                    difference,
                    allowed
                )
            )
        }
    }

    list(equal = !length(errors), errors = as.list(errors))
}

omicsParityValidateResult <- function(path) {
    result <- omicsParityReadJson(path, .OMICS_PARITY_RESULT_SCHEMA, "baseline result")
    required <- c(
        "run_id", "backend", "environment", "thread_controls", "sampling", "runs",
        "diagnostics", "determinism", "scenario_summaries", "summary"
    )
    missing <- setdiff(required, names(result))
    if (length(missing)) {
        stop(
            sprintf("Baseline result is missing: %s", paste(missing, collapse = ", ")),
            call. = FALSE
        )
    }
    if (!identical(result$backend, "memory")) {
        stop("Baseline result backend must be memory", call. = FALSE)
    }
    result
}

omicsParityQuantile <- function(values, probability) {
    values <- as.numeric(values)
    values <- values[is.finite(values)]
    if (!length(values)) {
        return(NA_real_)
    }
    unname(stats::quantile(values, probs = probability, names = FALSE, type = 8))
}

omicsParitySummarizeMeasurements <- function(runs) {
    successful <- Filter(\(run) identical(run$status, "passed"), runs)
    metric <- function(field) {
        values <- vapply(successful, \(run) {
            value <- run$metrics[[field]]
            if (is.null(value) || length(value) != 1L) NA_real_ else as.numeric(value)
        }, numeric(1))
        values <- values[is.finite(values)]
        list(
            median = if (length(values)) stats::median(values) else NA_real_,
            p95 = omicsParityQuantile(values, 0.95),
            maximum = if (length(values)) max(values) else NA_real_
        )
    }
    list(
        completed = length(successful),
        failed = length(runs) - length(successful),
        elapsed_seconds = metric("elapsed_seconds"),
        peak_tree_rss_bytes = metric("peak_tree_rss_bytes"),
        retained_tree_rss_bytes = metric("retained_tree_rss_bytes"),
        peak_disk_bytes = metric("peak_disk_bytes"),
        peak_artifact_disk_bytes = metric("peak_artifact_disk_bytes"),
        committed_input_bytes = metric("committed_input_bytes"),
        bounded_query_p95_seconds = metric("bounded_query_p95_seconds"),
        final_file_count = metric("final_file_count")
    )
}
