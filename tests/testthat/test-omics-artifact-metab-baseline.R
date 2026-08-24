metabBaseline029RepoPath <- function(...) {
    file.path(normalizePath(testthat::test_path("..", ".."), mustWork = TRUE), ...)
}

metabBaseline029ReadJson <- function(...) {
    jsonlite::read_json(metabBaseline029RepoPath(...), simplifyVector = FALSE)
}

metabBaseline029Named <- function(values, field) {
    stats::setNames(
        values,
        vapply(values, `[[`, character(1), field)
    )
}

test_that("metabolomics format support is exact and fail closed", {
    manifest <- metabBaseline029ReadJson(
        "tests", "testdata", "omics-parity", "metabolomics", "manifest.json"
    )
    decisions <- metabBaseline029Named(manifest$support_decisions, "input_format")
    catalogue <- Filter(
        \(entry) identical(entry$omic_type, "metabolomics"),
        workflowFormatSupportCatalogue()
    )
    catalogue <- metabBaseline029Named(catalogue, "ui_value")
    expect_setequal(names(decisions), names(catalogue))
    for (format in names(decisions)) {
        expect_identical(
            decisions[[format]]$support_status,
            catalogue[[format]]$support_status
        )
    }
    for (format in c("progenesis", "xcms", "compound_discoverer")) {
        expect_error(
            resolveWorkflowFormatSupport("metabolomics", format),
            class = "multischolar_format_detection_only"
        )
    }
    expect_identical(
        resolveWorkflowFormatSupport("metabolomics", "custom")$ui_value,
        "custom"
    )
    expect_identical(
        resolveWorkflowFormatSupport(
            "metabolomics",
            "msdial",
            detectedFormat = "msdial",
            detectionConfidence = 1
        )$ui_value,
        "msdial"
    )
    expect_error(
        resolveWorkflowFormatSupport("metabolomics", "unknown_vendor"),
        class = "multischolar_format_unknown"
    )
})

test_that("reviewed MS-DIAL and custom fixtures match independent oracles", {
    oracle <- metabBaseline029ReadJson(
        "tests", "testdata", "omics-parity", "metabolomics", "memory-oracle.json"
    )
    oracles <- metabBaseline029Named(oracle$oracles, "oracle_id")
    msdial <- suppressMessages(importMetabMSDIALData(metabBaseline029RepoPath(
        "tests", "testdata", "e2e", "metab_lc", "seed_lcms_pos.csv"
    )))
    expected <- oracles$msdial_lcms_pos_v1
    expect_identical(nrow(msdial$data), expected$rows)
    expect_identical(msdial$sample_columns, unlist(expected$sample_columns))
    expect_equal(
        sum(as.matrix(msdial$data[, msdial$sample_columns])),
        expected$quantity_sum,
        tolerance = 1e-8
    )
    expect_identical(as.character(msdial$data[["Alignment ID"]][[1L]]),
        expected$first_feature
    )
    expect_identical(msdial$data$Name[[1L]], expected$first_annotation)
    expect_identical(msdial$format, "msdial")

    manifest <- metabBaseline029ReadJson(
        "tests", "testdata", "omics-parity", "metabolomics", "manifest.json"
    )
    scenarios <- manifest$fixture_scenarios
    scenarios <- Filter(\(scenario) identical(scenario$format, "custom"), scenarios)
    for (scenario in scenarios) {
        expected <- oracles[[scenario$oracle_id]]$assays
        expected <- metabBaseline029Named(expected, "assay")
        paths <- unlist(scenario$data_paths, use.names = FALSE)
        for (index in seq_along(paths)) {
            data <- readr::read_tsv(
                metabBaseline029RepoPath(paths[[index]]),
                show_col_types = FALSE
            )
            assay <- scenario$assays[[index]]
            current <- expected[[assay]]
            quantities <- as.matrix(data[, 5:ncol(data)])
            expect_identical(nrow(data), current$rows)
            expect_identical(ncol(quantities), current$samples)
            expect_equal(sum(quantities), current$quantity_sum, tolerance = 1e-8)
            expect_identical(data[[1L]][[1L]], current$first_feature)
        }
    }
})

test_that("metabolomics fixture catalogue covers required edge identities", {
    fixture_dir <- metabBaseline029RepoPath(
        "tests", "testdata", "mci", "fx", "metabolomics"
    )
    files <- basename(list.files(fixture_dir, full.names = TRUE))
    required <- c(
        "happy-path", "bad-names", "duplicates", "missingness",
        "multi-assay-mixed", "nonfinite", "invalid-design",
        "large-enough-for-plots", "small-n", "multi-factor"
    )
    for (case in required) {
        expect_true(any(grepl(case, files, fixed = TRUE)), info = case)
    }
    nonfinite <- readr::read_tsv(
        file.path(fixture_dir, "metabolomics-lc-nonfinite_data.tsv"),
        show_col_types = FALSE
    )
    expect_true(any(!is.finite(as.matrix(nonfinite[, 6:ncol(nonfinite)]))))
    duplicates <- readr::read_tsv(
        file.path(fixture_dir, "metabolomics-lc-duplicates_data.tsv"),
        show_col_types = FALSE
    )
    expect_true(anyDuplicated(duplicates[[1L]]) > 0L)
    design <- readr::read_tsv(
        metabBaseline029RepoPath(
            "tests", "testdata", "e2e", "metab_combined", "design.tsv"
        ),
        show_col_types = FALSE
    )
    expect_identical(
        design$sample,
        c("WT_1", "WT_2", "WT_3", "KO_1", "KO_2", "KO_3")
    )
    expect_identical(unique(design$group), c("WT", "KO"))
    import_tests <- paste(readLines(metabBaseline029RepoPath(
        "tests", "testthat", "test-module-ci-metab-import.R"
    ), warn = FALSE), collapse = "\n")
    design_tests <- paste(readLines(metabBaseline029RepoPath(
        "tests", "testthat", "test-module-ci-metab-design.R"
    ), warn = FALSE), collapse = "\n")
    qc_tests <- paste(readLines(metabBaseline029RepoPath(
        "tests", "testthat", "test-module-ci-metab-qc.R"
    ), warn = FALSE), collapse = "\n")
    expect_match(import_tests, "non_numeric_samples", fixed = TRUE)
    expect_match(import_tests, "zero_missing", fixed = TRUE)
    expect_match(design_tests, "reordered assays", fixed = TRUE)
    expect_match(design_tests, "corrupt sample sets", fixed = TRUE)
    expect_match(qc_tests, "zeros, missingness", fixed = TRUE)
})

test_that("baseline freezes branch and workload claim boundaries", {
    manifest <- metabBaseline029ReadJson(
        "tests", "testdata", "omics-parity", "metabolomics", "manifest.json"
    )
    expect_setequal(names(manifest$branch_oracles), c(
        "qc", "normalization", "ruv", "da", "summary"
    ))
    expect_true(all(vapply(manifest$branch_oracles, length, integer(1)) > 0L))
    contracts <- unlist(manifest$workload_contracts, use.names = FALSE)
    expect_length(contracts, 5L)
    workload_ids <- vapply(contracts, \(relative) {
        jsonlite::read_json(
            metabBaseline029RepoPath(
                "tests", "testdata", "omics-parity", "metabolomics", relative
            ),
            simplifyVector = FALSE
        )$workload_id
    }, character(1))
    expect_identical(anyDuplicated(workload_ids), 0L)
    expect_true(any(grepl("public-ci", workload_ids, fixed = TRUE)))
    expect_true(any(grepl("public-representative", workload_ids, fixed = TRUE)))
    expect_match(
        manifest$claim_boundary,
        "vendor support and biological parity require reviewed fixtures",
        fixed = TRUE
    )
})

test_that("representative memory resource evidence has frozen finite gates", {
    path <- metabBaseline029RepoPath(
        "tests", "testdata", "omics-parity", "metabolomics", "resources",
        "mixed-memory-baseline-v1.json"
    )
    testthat::skip_if_not(file.exists(path), "resource evidence not generated")
    evidence <- jsonlite::read_json(path, simplifyVector = FALSE)
    expect_identical(evidence$status, "passed")
    expect_identical(evidence$summary$completed, 5L)
    expect_length(evidence$runs, 5L)
    expect_true(all(unlist(evidence$determinism, use.names = FALSE)))
    gates <- evidence$release_gates
    expect_true(all(vapply(gates, \(value) {
        is.numeric(value) && length(value) == 1L && is.finite(value) && value >= 0
    }, logical(1))))
    for (run in evidence$runs) {
        expect_identical(run$status, "passed")
        expect_true(run$worker$workflow_evidence$truth_valid)
        expect_identical(run$worker$session_evidence$resource_count, 0L)
        expect_identical(run$worker$report_evidence$file_count, 0L)
        expect_false(run$worker$native_resources$arrow$loaded)
        expect_false(run$worker$native_resources$duckdb$loaded)
    }
})
