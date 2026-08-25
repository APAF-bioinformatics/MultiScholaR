lipidBaselineRepoPath <- function(...) {
    file.path(normalizePath(testthat::test_path("..", ".."), mustWork = TRUE), ...)
}

lipidBaselineManifest <- function() {
    jsonlite::read_json(
        lipidBaselineRepoPath(
            "tests", "testdata", "omics-parity", "lipidomics", "manifest.json"
        ),
        simplifyVector = TRUE
    )
}

lipidBaselineOracle <- function() {
    jsonlite::read_json(
        lipidBaselineRepoPath(
            "tests", "testdata", "omics-parity", "lipidomics",
            "memory-oracle.json"
        ),
        simplifyVector = TRUE
    )
}

test_that("lipidomics support decisions fail closed", {
    manifest <- lipidBaselineManifest()
    decisions <- manifest$support_decisions
    expect_identical(
        decisions$support_status[decisions$input_format == "lipidsearch"],
        "scientifically_supported"
    )
    expect_identical(
        decisions$support_status[decisions$input_format == "msdial"],
        "reader_characterized"
    )
    expect_identical(
        decisions$support_status[decisions$input_format == "custom"],
        "reader_characterized"
    )
    blocked <- c("progenesis", "xcms", "compound_discoverer")
    expect_true(all(
        decisions$support_status[decisions$input_format %in% blocked] ==
            "detection_only"
    ))
    for (format in blocked) {
        expect_error(
            resolveWorkflowFormatSupport(
                "lipidomics",
                format,
                format,
                0.9
            ),
            class = "multischolar_format_detection_only"
        )
    }
})

test_that("LipidSearch fixtures match independent identity and arithmetic oracles", {
    manifest <- lipidBaselineManifest()
    oracle <- lipidBaselineOracle()
    scenarios <- manifest$fixture_scenarios[1:3, ]
    for (index in seq_len(nrow(scenarios))) {
        scenario <- scenarios[index, ]
        expected <- oracle$oracles[oracle$oracles$oracle_id == scenario$oracle_id, ]
        contract <- jsonlite::read_json(
            lipidBaselineRepoPath(
                "tests", "testdata", "omics-parity", "lipidomics",
                manifest$workload_contracts[[index]]
            ),
            simplifyVector = TRUE
        )
        contract_oracle <- contract$adapter_parameters$oracle
        path <- lipidBaselineRepoPath(scenario$data_paths[[1L]][[1L]])
        imported <- suppressMessages(importLipidSearchData(path))
        data <- imported$data
        sample_columns <- imported$sample_columns
        expect_identical(nrow(data), expected$rows)
        expect_identical(sample_columns, unlist(expected$sample_columns[[1L]]))
        expect_equal(sum(as.matrix(data[sample_columns])), expected$quantity_sum)
        expect_identical(data$LipidName[[1L]], expected$first_lipid)
        expect_identical(data$LipidClass[[1L]], expected$first_class)
        expect_identical(data$IonType[[1L]], expected$first_adduct)
        expect_identical(contract_oracle$rows, expected$rows)
        expect_identical(
            contract_oracle$sample_columns,
            unlist(expected$sample_columns[[1L]])
        )
        expect_identical(contract_oracle$quantity_sum, expected$quantity_sum)
        expect_identical(contract_oracle$first_lipid, expected$first_lipid)
        expect_identical(contract_oracle$first_class, expected$first_class)
        expect_identical(contract_oracle$first_adduct, expected$first_adduct)
    }
})

test_that("mixed LipidSearch S4 state preserves lipid and assay identity", {
    files <- c(
        LCMS_Pos = "lipidsearch_lcms_pos.txt",
        LCMS_Neg = "lipidsearch_lcms_neg.txt",
        GCMS = "lipidsearch_gcms.txt"
    )
    assays <- lapply(files, \(file) suppressMessages(importLipidSearchData(
        lipidBaselineRepoPath("tests", "testdata", "e2e", "lipid_canonical", file)
    ))$data)
    samples <- grep("^(WT|KO)_", names(assays[[1L]]), value = TRUE)
    groups <- ifelse(grepl("^KO", samples), "KO", "WT")
    design <- data.frame(
        Run = samples,
        group = groups,
        tech_rep_group = paste0(groups, "_", ave(seq_along(groups), groups)),
        stringsAsFactors = FALSE
    )
    object <- createLipidomicsAssayData(
        assays,
        design,
        lipid_id_column = "LipidName",
        annotation_id_column = "LipidClass",
        sample_id = "Run",
        group_id = "group",
        technical_replicate_id = "tech_rep_group",
        internal_standard_regex = "^IS_"
    )
    expect_s4_class(object, "LipidomicsAssayData")
    expect_identical(names(object@lipid_data), names(files))
    expect_identical(methods::validObject(object), TRUE)
    expect_identical(object@lipid_data$LCMS_Pos$LipidName[[1L]], "PC 34:1")
    expect_identical(object@lipid_data$LCMS_Neg$LipidClass[[1L]], "PI")
    expect_identical(object@lipid_data$GCMS$IonType[[1L]], "[M-H]-")
})

test_that("baseline declares every downstream scientific branch", {
    branches <- lipidBaselineManifest()$branch_oracles
    expect_setequal(names(branches), c(
        "qc", "normalization", "ruv", "da", "summary"
    ))
    expect_true(all(vapply(branches, length, integer(1)) > 0L))
    expect_match(
        lipidBaselineManifest()$claim_boundary,
        "vendor support and biological parity require reviewed fixtures",
        fixed = TRUE
    )
})

test_that("fixture and representative resource gates are finite and frozen", {
    manifest <- lipidBaselineManifest()
    evidence <- manifest$resource_evidence
    expect_identical(nrow(evidence), 7L)
    workload_ids <- character()
    for (index in seq_len(nrow(evidence))) {
        entry <- evidence[index, ]
        contract <- jsonlite::read_json(
            lipidBaselineRepoPath(
                "tests", "testdata", "omics-parity", "lipidomics",
                entry$contract
            ),
            simplifyVector = FALSE
        )
        resource <- jsonlite::read_json(
            lipidBaselineRepoPath(
                "tests", "testdata", "omics-parity", "lipidomics",
                entry$result
            ),
            simplifyVector = FALSE
        )
        expect_identical(resource$status, "passed", info = entry$result)
        expect_identical(resource$summary$completed, 5L, info = entry$result)
        expect_identical(
            resource$workload_id,
            contract$workload_id,
            info = entry$result
        )
        expect_identical(
            resource$preparation$metadata$evidence_class,
            entry$evidence_class,
            info = entry$result
        )
        expect_true(all(vapply(resource$runs, \(run) {
            identical(
                run$worker$workflow_evidence$evidence_class,
                entry$evidence_class
            )
        }, logical(1))), info = entry$result)
        gates <- unlist(resource$release_gates, use.names = TRUE)
        expect_setequal(names(gates), c(
            "maximum_peak_tree_rss_bytes",
            "maximum_retained_tree_rss_bytes",
            "maximum_elapsed_seconds",
            "maximum_peak_disk_bytes",
            "maximum_transient_disk_bytes",
            "maximum_peak_artifact_disk_bytes",
            "maximum_peak_staging_snapshot_disk_bytes",
            "maximum_peak_duckdb_spill_bytes",
            "maximum_peak_committed_disk_bytes",
            "maximum_peak_final_disk_bytes",
            "maximum_final_disk_bytes",
            "maximum_final_file_count",
            "maximum_committed_input_bytes",
            "maximum_committed_input_file_count",
            "maximum_session_resource_count",
            "maximum_report_file_count",
            "maximum_bounded_query_p95_seconds"
        ))
        expect_true(all(is.finite(gates)), info = entry$result)
        expect_true(all(gates >= 0), info = entry$result)
        positive <- c(
            "maximum_peak_tree_rss_bytes",
            "maximum_retained_tree_rss_bytes",
            "maximum_elapsed_seconds",
            "maximum_peak_disk_bytes",
            "maximum_peak_final_disk_bytes",
            "maximum_final_disk_bytes",
            "maximum_final_file_count",
            "maximum_committed_input_bytes",
            "maximum_committed_input_file_count",
            "maximum_bounded_query_p95_seconds"
        )
        expect_true(all(gates[positive] > 0), info = entry$result)
        observed_zero <- c(
            "maximum_peak_artifact_disk_bytes",
            "maximum_peak_staging_snapshot_disk_bytes",
            "maximum_peak_duckdb_spill_bytes",
            "maximum_peak_committed_disk_bytes",
            "maximum_session_resource_count",
            "maximum_report_file_count"
        )
        expect_identical(
            unname(gates[observed_zero]),
            rep(0, length(observed_zero)),
            info = entry$result
        )
        expect_match(
            resource$gate_policy$claim_boundary,
            "resource",
            fixed = TRUE,
            info = entry$result
        )
        workload_ids <- c(workload_ids, resource$workload_id)
    }
    expect_identical(anyDuplicated(workload_ids), 0L)
})
