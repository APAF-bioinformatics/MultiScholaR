.protNondiaRepoPath <- function(...) {
    file.path(normalizePath(testthat::test_path("..", ".."), mustWork = TRUE), ...)
}

.protNondiaReadJson <- function(name) {
    jsonlite::read_json(
        .protNondiaRepoPath(
            "tests",
            "testdata",
            "omics-parity",
            "proteomics",
            name
        ),
        simplifyVector = FALSE
    )
}

.protNondiaNamed <- function(values, field) {
    keys <- vapply(values, \(value) value[[field]], character(1))
    stats::setNames(values, keys)
}

.protNondiaImporter <- function(format) {
    switch(
        format,
        maxquant = importMaxQuantData,
        fragpipe = importFragPipeData,
        pd_tmt = importProteomeDiscovererTMTData
    )
}

.protNondiaCanonicalDigest <- function(data, keys) {
    data <- as.data.frame(data)
    ordering <- do.call(order, c(data[keys], list(method = "radix")))
    data <- data[ordering, , drop = FALSE]
    data[] <- lapply(data, \(column) {
        if (is.numeric(column)) {
            return(sprintf("%.10f", column))
        }
        as.character(column)
    })
    payload <- jsonlite::toJSON(
        data,
        dataframe = "rows",
        na = "null",
        digits = 17
    )
    digest::digest(payload, algo = "sha256", serialize = FALSE)
}

.protNondiaImport <- function(scenario) {
    importer <- .protNondiaImporter(scenario$format)
    suppressMessages(importer(.protNondiaRepoPath(scenario$fixture_path)))
}

.protNondiaBuildObject <- function(imported, workflow_type) {
    data <- as.data.frame(imported$data)
    mapping <- imported$column_mapping
    runs <- unique(as.character(data[[mapping$run_col]]))
    groups <- ifelse(grepl("KO", runs, fixed = TRUE), "KO", "WT")
    replicates <- ave(seq_along(runs), groups, FUN = seq_along)
    design <- data.frame(
        Run = runs,
        group = groups,
        replicates = paste0("R", replicates),
        stringsAsFactors = FALSE
    )
    manager <- new.env(parent = emptyenv())
    manager$saveState <- \(state_name, s4_data_object, ...) {
        manager$object <- s4_data_object
        invisible(state_name)
    }
    workflow <- list2env(list(
        data_cln = data,
        column_mapping = mapping,
        design_matrix = design,
        config_list = list(
            globalParameters = list(workflow_type = workflow_type)
        ),
        state_manager = manager
    ), parent = emptyenv())
    suppressMessages(buildProtDesignStateCheckpoint(
        workflowData = workflow,
        workflowType = workflow_type,
        actionLabel = "non-DIA baseline",
        validateColumnMapping = TRUE
    ))
    manager$object
}

.protNondiaNormalize <- function(object) {
    invisible(utils::capture.output(
        normalized <- suppressMessages(normaliseBetweenSamples(
            object,
            normalisation_method = "cyclicloess"
        ))
    ))
    normalized
}

.protNondiaRunDa <- function(object) {
    matrix <- as.matrix(object@protein_quant_table[
        ,
        setdiff(names(object@protein_quant_table), object@protein_id_column),
        drop = FALSE
    ])
    rownames(matrix) <- object@protein_quant_table[[object@protein_id_column]]
    design <- object@design_matrix
    rownames(design) <- design$Run
    results <- suppressMessages(runTestsContrasts(
        data = matrix,
        contrast_strings = "KO_vs_WT=groupKO-groupWT",
        design_matrix = design,
        formula_string = "~ 0 + group",
        treat_lfc_cutoff = 0,
        eBayes_trend = TRUE,
        eBayes_robust = TRUE
    ))$results[[1L]]
    results$feature_id <- rownames(results)
    results[order(results$feature_id, method = "radix"), , drop = FALSE]
}

.protNondiaRawLog2Fc <- function(imported) {
    data <- as.data.frame(imported$data)
    mapping <- imported$column_mapping
    features <- as.character(data[[mapping$protein_col]])
    samples <- as.character(data[[mapping$run_col]])
    values <- log2(as.numeric(data[[mapping$quantity_col]]) + 1)
    means <- function(selected) {
        ids <- features[selected]
        sums <- rowsum(values[selected], ids, reorder = FALSE)
        counts <- rowsum(rep.int(1, length(ids)), ids, reorder = FALSE)
        stats::setNames(as.numeric(sums / counts), rownames(sums))
    }
    means(grepl("KO", samples, fixed = TRUE)) -
        means(grepl("WT", samples, fixed = TRUE))
}

test_that("non-DIA support and reachable-stage decisions are explicit", {
    manifest <- .protNondiaReadJson("manifest.json")
    decisions <- .protNondiaNamed(manifest$support_decisions, "capability_id")
    expect_identical(manifest$backend, "memory")
    expect_setequal(names(decisions), c(
        "proteomics.maxquant.protein.lfq.v1",
        "proteomics.fragpipe.protein.lfq.v1",
        "proteomics.pd_tmt.protein.tmt.v1",
        "proteomics.spectronaut.protein.lfq.v1",
        "proteomics.spectronaut.peptide.lfq.v1"
    ))

    supported <- decisions[vapply(
        decisions,
        \(decision) identical(decision$support_status, "scientifically_supported"),
        logical(1)
    )]
    for (decision in supported) {
        stages <- unlist(decision$reachable_stages, use.names = FALSE)
        expect_true("protein_qc_bypass" %in% stages)
        expect_false(any(grepl("^protein_qc$", stages)))
        expect_true(length(decision$e2e_lane) == 1L)
    }

    blocked <- decisions[vapply(
        decisions,
        \(decision) identical(decision$support_status, "advertised_unverified"),
        logical(1)
    )]
    expect_length(blocked, 2L)
    expect_true(all(vapply(blocked, \(decision) {
        is.null(decision$e2e_lane) && !length(decision$reachable_stages)
    }, logical(1))))
})

test_that("reviewed vendor fixtures reproduce exact import oracles", {
    manifest <- .protNondiaReadJson("manifest.json")
    oracle <- .protNondiaReadJson("memory-oracle.json")
    expected <- .protNondiaNamed(oracle$imports, "scenario_id")
    tolerance <- manifest$comparison_contract$import_numeric_tolerance

    for (scenario in manifest$fixture_scenarios) {
        imported <- .protNondiaImport(scenario)
        data <- as.data.frame(imported$data)
        mapping <- imported$column_mapping
        observed <- list(
            data_type = imported$data_type,
            mapping = mapping,
            rows = nrow(data),
            columns = ncol(data),
            column_names = as.list(names(data)),
            feature_ids = as.list(unique(as.character(
                data[[mapping$protein_col]]
            ))),
            sample_ids = as.list(unique(as.character(
                data[[mapping$run_col]]
            ))),
            table_sha256 = .protNondiaCanonicalDigest(
                data,
                c(mapping$protein_col, mapping$run_col)
            )
        )
        current <- expected[[scenario$scenario_id]]
        exact_fields <- unlist(
            manifest$comparison_contract$import_exact_fields,
            use.names = FALSE
        )
        for (field in exact_fields) {
            expect_identical(observed[[field]], current[[field]], info = paste(
                scenario$scenario_id,
                field
            ))
        }
        quantity <- as.numeric(data[[mapping$quantity_col]])
        expect_equal(
            sum(quantity),
            current$quantity_sum,
            tolerance = tolerance$absolute
        )
        expect_equal(
            min(quantity),
            current$quantity_min,
            tolerance = tolerance$absolute
        )
        expect_equal(
            max(quantity),
            current$quantity_max,
            tolerance = tolerance$absolute
        )
    }
})

test_that("current design, normalization, and DA reproduce the scientific oracle", {
    manifest <- .protNondiaReadJson("manifest.json")
    oracle <- .protNondiaReadJson("memory-oracle.json")
    scientific <- .protNondiaNamed(
        oracle$scientific_oracles,
        "scientific_oracle_id"
    )[["shared_5x6_cyclicloess_v1"]]
    expected_fc <- unlist(scientific$raw_log2fc, use.names = TRUE)
    expected_normalized <- lapply(
        scientific$normalized_values,
        \(values) as.numeric(unlist(values, use.names = FALSE))
    )
    expected_da <- .protNondiaNamed(scientific$da_rows, "feature_id")
    norm_tolerance <- manifest$comparison_contract$normalization_numeric_tolerance
    da_tolerance <- manifest$comparison_contract$da_numeric_tolerance

    for (scenario in manifest$fixture_scenarios) {
        imported <- .protNondiaImport(scenario)
        raw_fc <- .protNondiaRawLog2Fc(imported)
        expect_equal(
            raw_fc[names(expected_fc)],
            expected_fc,
            tolerance = norm_tolerance$absolute,
            info = scenario$scenario_id
        )
        object <- .protNondiaBuildObject(
            imported,
            if (identical(scenario$format, "pd_tmt")) "TMT" else "LFQ"
        )
        normalized <- .protNondiaNormalize(object)
        expected_columns <- c(
            "Protein.Ids",
            unique(as.character(
                imported$data[[imported$column_mapping$run_col]]
            ))
        )
        expect_identical(names(normalized@protein_quant_table), expected_columns)
        for (feature_id in names(expected_normalized)) {
            observed <- normalized@protein_quant_table[
                normalized@protein_quant_table$Protein.Ids == feature_id,
                -1L,
                drop = TRUE
            ]
            expect_equal(
                as.numeric(observed),
                expected_normalized[[feature_id]],
                tolerance = norm_tolerance$absolute,
                info = paste(scenario$scenario_id, feature_id)
            )
        }

        da <- .protNondiaRunDa(normalized)
        expect_identical(
            names(da),
            unlist(scientific$da_column_names, use.names = FALSE)
        )
        for (feature_id in names(expected_da)) {
            observed <- da[da$feature_id == feature_id, , drop = FALSE]
            expected <- expected_da[[feature_id]]
            expect_identical(observed$feature_id, feature_id)
            expect_true(is.na(observed$fdr_qvalue))
            for (field in c(
                "logFC",
                "AveExpr",
                "t",
                "raw_pvalue",
                "adj.P.Val",
                "fdr_value_bh_adjustment"
            )) {
                expect_equal(
                    observed[[field]],
                    expected[[field]],
                    tolerance = da_tolerance$absolute,
                    info = paste(scenario$scenario_id, feature_id, field)
                )
            }
        }
        expect_identical(
            da$feature_id[da$fdr_value_bh_adjustment <= 0.05],
            character()
        )
    }
})

test_that("E2E metadata records bypass, report, and normalization branches", {
    e2e <- jsonlite::read_json(
        .protNondiaRepoPath("tests", "testdata", "e2e", "manifest.json"),
        simplifyVector = FALSE
    )
    lanes <- .protNondiaNamed(e2e$lanes, "lane_id")
    expected <- c(
        prot_lfq = "maxquant_lfq",
        prot_lfq_fragpipe = "fragpipe_lfq",
        prot_tmt = "pd_tmt"
    )
    for (lane_id in names(expected)) {
        lane <- lanes[[lane_id]]
        expect_identical(lane$scientific_baseline_id, expected[[lane_id]])
        expect_identical(lane$qc_route, "top_level_protein_qc_bypass")
        expect_identical(lane$normalization_method, "cyclicloess")
        expect_true("qc_bypass" %in% unlist(lane$touchpoints, use.names = FALSE))
        expect_false(any(grepl(
            "^qc_",
            unlist(lane$module_families, use.names = FALSE)
        )))
        expect_true(lane$report_template %in% c(
            "LFQ_report.rmd",
            "TMT_report.rmd"
        ))
    }
})

test_that("tuple-specific resource evidence meets frozen baseline requirements", {
    manifest <- .protNondiaReadJson("manifest.json")
    margin <- manifest$resource_gate_policy$materiality_margin_fraction
    resource_dir <- .protNondiaRepoPath(
        "tests",
        "testdata",
        "omics-parity",
        "proteomics",
        "resources"
    )
    paths <- sort(list.files(
        resource_dir,
        pattern = "[.]json$",
        full.names = TRUE
    ))
    expect_length(paths, 3L)
    capability_ids <- character()
    for (path in paths) {
        evidence <- jsonlite::read_json(path, simplifyVector = FALSE)
        capability_ids <- c(
            capability_ids,
            evidence$capability$capability_id
        )
        expect_identical(evidence$status, "passed", info = basename(path))
        expect_identical(
            evidence$release_gates$materiality_margin_fraction,
            margin
        )
        expect_gte(
            evidence$summary$completed,
            manifest$resource_gate_policy$minimum_repetitions_per_tuple
        )
        expect_true(all(unlist(evidence$determinism, use.names = FALSE)))
        expect_identical(evidence$preparation$fresh_processes, 2L)
        expect_true(evidence$preparation$deterministic)
        expect_length(evidence$runs, 5L)

        gate_metrics <- c(
            maximum_peak_tree_rss_bytes = "peak_tree_rss_bytes",
            maximum_retained_tree_rss_bytes = "retained_tree_rss_bytes",
            maximum_elapsed_seconds = "elapsed_seconds",
            maximum_peak_disk_bytes = "peak_disk_bytes",
            maximum_final_disk_bytes = "final_disk_bytes",
            maximum_final_file_count = "final_file_count"
        )
        for (gate_name in names(gate_metrics)) {
            metric_name <- gate_metrics[[gate_name]]
            observed <- vapply(
                evidence$runs,
                \(run) run$metrics[[metric_name]],
                numeric(1)
            )
            expected_gate <- max(observed) * (1 + margin)
            if (grepl("seconds$", gate_name)) {
                expected_gate <- ceiling(expected_gate * 1000) / 1000
            } else {
                expected_gate <- ceiling(expected_gate)
            }
            expect_equal(
                evidence$release_gates[[gate_name]],
                expected_gate,
                tolerance = 1e-12,
                info = paste(basename(path), gate_name)
            )
        }

        query_p95 <- vapply(
            evidence$runs,
            \(run) run$worker$query_evidence[[1L]]$p95_seconds,
            numeric(1)
        )
        expected_query_gate <- ceiling(
            max(query_p95) * (1 + margin) * 1000
        ) / 1000
        expect_equal(
            evidence$release_gates$maximum_bounded_query_p95_seconds,
            expected_query_gate,
            tolerance = 1e-12
        )
        for (gate_name in c(
            "maximum_peak_artifact_disk_bytes",
            "maximum_duckdb_spill_bytes",
            "maximum_session_resource_count",
            "maximum_report_file_count"
        )) {
            expect_identical(evidence$release_gates[[gate_name]], 0L)
        }

        for (run in evidence$runs) {
            expect_identical(run$status, "passed")
            expect_true(run$retention_marker_observed)
            for (gate_name in names(gate_metrics)) {
                metric_name <- gate_metrics[[gate_name]]
                expect_true(
                    run$metrics[[metric_name]] <=
                        evidence$release_gates[[gate_name]],
                    info = paste(basename(path), gate_name)
                )
            }
            expect_gt(run$metrics$committed_input_bytes, 0)
            expect_identical(run$metrics$committed_input_file_count, 2L)
            expect_identical(run$metrics$peak_artifact_disk_bytes, 0L)
            expect_identical(
                run$metrics$peak_disk_category_bytes$duckdb_spill,
                0L
            )
            expect_lte(
                run$worker$query_evidence[[1L]]$p95_seconds,
                evidence$release_gates$maximum_bounded_query_p95_seconds
            )
            expect_true(run$worker$workflow_evidence$truth_valid)
            expect_true(
                run$worker$workflow_evidence$differential_direction$valid
            )
            expect_identical(run$worker$session_evidence$resource_count, 0L)
            expect_identical(run$worker$report_evidence$file_count, 0L)
            expect_false(run$worker$native_resources$arrow$loaded)
            expect_false(run$worker$native_resources$duckdb$loaded)
            expect_identical(
                run$binding$contract_sha256,
                evidence$contract_sha256
            )
            expect_true(run$validation$worker$valid)
            expect_true(run$validation$finite_release_metrics)
        }

        expect_length(evidence$diagnostics, 1L)
        diagnostic <- evidence$diagnostics[[1L]]
        expect_identical(
            diagnostic$measurement_class,
            "allocation_copy_diagnostic"
        )
        expect_true(diagnostic$worker$allocation_diagnostics$available)
        expect_gt(
            diagnostic$worker$allocation_diagnostics$allocated_bytes,
            0
        )
        expect_gt(
            diagnostic$worker$allocation_diagnostics$allocation_records,
            0L
        )
        expect_match(
            diagnostic$worker$allocation_diagnostics$profile_sha256,
            "^[[:xdigit:]]{64}$"
        )
    }
    expect_setequal(capability_ids, c(
        "proteomics.maxquant.protein.lfq.v1",
        "proteomics.fragpipe.protein.lfq.v1",
        "proteomics.pd_tmt.protein.tmt.v1"
    ))
})
