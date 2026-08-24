.metabCodecRepoPath <- function(...) {
    file.path(normalizePath(testthat::test_path("..", ".."), mustWork = TRUE), ...)
}

.metabCodecLogicalKey <- function(
    workflow_slug = "custom",
    generation_id = "generation-001",
    stage_id = "design",
    state_role = "assay_data"
) {
    list(
        project_id = "metabolomics-codec-project",
        omic_type = "metabolomics",
        workflow_slug = workflow_slug,
        stage_id = stage_id,
        state_role = state_role,
        generation_id = generation_id
    )
}

.metabCodecWriteReadPayloads <- function(bundle) {
    payloads <- Map(
        \(payload, metadata, index) {
            encoded <- structure(
                list(payload = payload, metadata = metadata),
                class = c("MultiScholaRArtifactRectangular", "list")
            )
            path <- tempfile(sprintf("metab-codec-%03d-", index), fileext = ".parquet")
            do.call(
                arrow::write_parquet,
                c(
                    list(x = payload, sink = path),
                    artifactParquetWriteArgs(encoded)
                )
            )
            arrow::read_parquet(path, as_data_frame = FALSE)
        },
        bundle$payloads,
        bundle$metadata$payloads,
        seq_along(bundle$payloads)
    )
    names(payloads) <- names(bundle$payloads)
    bundle$payloads <- payloads
    bundle
}

.metabCodecTable <- function(prefix, offset = 0) {
    data.frame(
        feature_id = paste0(prefix, c("_003", "_001", "_002", "_004")),
        annotation = c("alpha/beta", "quoted 'name'", NA, "delta"),
        internal_standard = c(FALSE, TRUE, FALSE, FALSE),
        mz = c(101.1, 202.2, NA_real_, 404.4),
        `sample A` = c(0, 10 + offset, NA_real_, NaN),
        `sample/B` = c(Inf, 20 + offset, 30 + offset, -Inf),
        `sample-03` = c(3 + offset, 4 + offset, 5 + offset, 6 + offset),
        sample_04 = c(7 + offset, 8 + offset, 9 + offset, 10 + offset),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
}

.metabCodecAssay <- function(assay_names) {
    design <- data.frame(
        sample_key = c("sample-03", "sample A", "sample_04", "sample/B"),
        group = factor(
            c("treated", "control", "treated", "control"),
            levels = c("control", "treated")
        ),
        batch = c("B2", "B1", "B2", "B1"),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    assays <- lapply(seq_along(assay_names), \(index) {
        .metabCodecTable(assay_names[[index]], offset = index * 100)
    })
    names(assays) <- assay_names
    createMetaboliteAssayData(
        metabolite_data = assays,
        design_matrix = design,
        metabolite_id_column = "feature_id",
        annotation_id_column = "annotation",
        sample_id = "sample_key",
        group_id = "group",
        technical_replicate_id = NA_character_,
        database_identifier_type = "Mixed_CHEBI_Unknown",
        internal_standard_regex = "(^|_)IS($|_)",
        args = list(
            formula_string = "~ 0 + group + batch",
            thresholds = c(q_value = 0.05, log2_fc = 1),
            normalization = list(method = "none", centered = FALSE),
            quality_tier = factor("reviewed", levels = c("reviewed", "unused")),
            reference_date = as.Date("2026-08-17"),
            contrast_matrix = matrix(
                c(1, -1),
                nrow = 2L,
                dimnames = list(c("treated", "control"), "treated-control")
            )
        )
    )
}

.metabCodecExpectExactSlots <- function(before, after) {
    expect_identical(class(after), class(before))
    for (slot_name in methods::slotNames(before)) {
        expect_identical(
            methods::slot(after, slot_name),
            methods::slot(before, slot_name),
            info = slot_name
        )
    }
    expect_identical(methods::validObject(after, test = TRUE), TRUE)
}

.metabCodecDaResult <- function(assay) {
    sample_columns <- c("sample A", "sample/B", "sample-03", "sample_04")
    quantities <- as.matrix(assay@metabolite_data[[1L]][, sample_columns])
    quantities[!is.finite(quantities)] <- NA_real_
    quantities[is.na(quantities)] <- 0
    design <- stats::model.matrix(
        ~ 0 + group,
        data = assay@design_matrix[match(sample_columns, assay@design_matrix$sample_key), ]
    )
    fit <- limma::lmFit(quantities, design)
    contrasts <- list(
        treated_vs_control = data.frame(
            feature_id = assay@metabolite_data[[1L]]$feature_id,
            comparison = "treated_vs_control",
            log2FC = c(1.5, -1.2, 0.1, NA_real_),
            raw_pvalue = c(0.001, 0.02, 0.8, NA_real_),
            fdr_qvalue = c(0.004, 0.04, 0.8, NA_real_),
            stringsAsFactors = FALSE
        ),
        control_vs_treated = data.frame(
            feature_id = assay@metabolite_data[[1L]]$feature_id,
            comparison = "control_vs_treated",
            log2FC = c(-1.5, 1.2, -0.1, NA_real_),
            raw_pvalue = c(0.001, 0.02, 0.8, NA_real_),
            fdr_qvalue = c(0.004, 0.04, 0.8, NA_real_),
            stringsAsFactors = FALSE
        )
    )
    long <- do.call(rbind, unname(contrasts))
    row.names(long) <- NULL
    plot <- ggplot2::ggplot(
        contrasts[[1L]],
        ggplot2::aes(x = log2FC, y = -log10(fdr_qvalue))
    ) + ggplot2::geom_point()
    methods::new(
        "MetabolomicsDifferentialAbundanceResults",
        theObject = assay,
        fit.eb = fit,
        contrasts_results_table = contrasts,
        num_sig_diff_exp_bar_plot = list(treated_vs_control = plot),
        num_sig_diff_table = data.frame(
            comparison = names(contrasts),
            significant = c(2L, 2L),
            stringsAsFactors = FALSE
        ),
        volcano_plot = list(treated_vs_control = plot),
        interactive_volcano_plot = list(
            treated_vs_control = structure(
                list(x = list(data = contrasts[[1L]])),
                class = c("plotly", "htmlwidget")
            )
        ),
        p_value_dist_plot = list(treated_vs_control = plot),
        results_table_long = long,
        results_table_wide = merge(
            contrasts[[1L]],
            contrasts[[2L]],
            by = "feature_id",
            sort = FALSE,
            suffixes = c(".forward", ".reverse")
        )
    )
}

.metabCodecExpectDaSlots <- function(before, after) {
    exact_slots <- c(
        "theObject", "fit.eb", "contrasts_results_table", "num_sig_diff_table",
        "results_table_long", "results_table_wide"
    )
    for (slot_name in exact_slots) {
        expect_identical(
            methods::slot(after, slot_name),
            methods::slot(before, slot_name),
            info = slot_name
        )
    }
    plot_slots <- setdiff(methods::slotNames(before), exact_slots)
    for (slot_name in plot_slots) {
        expect_true(isTRUE(all.equal(
            methods::slot(after, slot_name),
            methods::slot(before, slot_name)
        )), info = slot_name)
    }
    expect_s4_class(after, "MetabolomicsDifferentialAbundanceResults")
    expect_identical(methods::validObject(after, test = TRUE), TRUE)
}

test_that("metabolomics assay codecs preserve LC, GC, mixed, and reordered state", {
    scenarios <- list(
        lc_only = "LCMS_Pos",
        gc_only = "GCMS",
        mixed = c("LCMS_Pos", "LCMS_Neg", "GCMS"),
        reordered = c("GCMS", "LCMS_Neg", "LCMS_Pos")
    )
    for (scenario in names(scenarios)) {
        before <- .metabCodecAssay(scenarios[[scenario]])
        bundle <- dehydrateMetabolomicsS4Artifact(
            before,
            .metabCodecLogicalKey(generation_id = paste0("generation-", scenario))
        )
        after <- hydrateMetabolomicsS4Artifact(
            .metabCodecWriteReadPayloads(bundle)
        )
        .metabCodecExpectExactSlots(before, after)
        expect_identical(names(after@metabolite_data), scenarios[[scenario]])
        expect_identical(
            vapply(
                bundle$metadata$identity_binding$assays,
                `[[`,
                character(1),
                "name"
            ),
            scenarios[[scenario]]
        )
        expect_identical(
            length(bundle$metadata$slot_values$metabolite_data$values),
            length(scenarios[[scenario]])
        )
    }
})

test_that("metabolomics assay identity is workflow and generation scoped", {
    assay <- .metabCodecAssay(c("LCMS_Pos", "GCMS"))
    first <- dehydrateMetabolomicsS4Artifact(
        assay,
        .metabCodecLogicalKey(workflow_slug = "custom", generation_id = "g1")
    )
    repeat_first <- dehydrateMetabolomicsS4Artifact(
        assay,
        .metabCodecLogicalKey(workflow_slug = "custom", generation_id = "g1")
    )
    other_generation <- dehydrateMetabolomicsS4Artifact(
        assay,
        .metabCodecLogicalKey(workflow_slug = "custom", generation_id = "g2")
    )
    other_workflow <- dehydrateMetabolomicsS4Artifact(
        assay,
        .metabCodecLogicalKey(workflow_slug = "msdial", generation_id = "g1")
    )
    ids <- \(bundle) vapply(
        bundle$metadata$identity_binding$assays,
        `[[`,
        character(1),
        "assay_id"
    )
    expect_identical(ids(first), ids(repeat_first))
    expect_false(any(ids(first) %in% ids(other_generation)))
    expect_false(any(ids(first) %in% ids(other_workflow)))
    expect_identical(anyDuplicated(ids(first)), 0L)
})

test_that("metabolomics DA codec composes exact tables and immutable R roles", {
    before <- .metabCodecDaResult(
        .metabCodecAssay(c("GCMS", "LCMS_Pos"))
    )
    key <- .metabCodecLogicalKey(
        generation_id = "da-generation-001",
        stage_id = "differential_abundance",
        state_role = "complete_results"
    )
    bundle <- dehydrateMetabolomicsS4Artifact(before, key)
    expect_s3_class(
        validateMetabolomicsS4Bundle(bundle),
        "MultiScholaRArtifactS4Bundle"
    )
    after <- hydrateMetabolomicsS4Artifact(
        .metabCodecWriteReadPayloads(bundle)
    )
    .metabCodecExpectDaSlots(before, after)
    nodes <- bundle$metadata$slot_values
    expect_identical(nodes$theObject$node_type, "s4_exact_role")
    for (slot_name in names(.artifactMetabSerializedSlots)) {
        expect_identical(nodes[[slot_name]]$node_type, "immutable_r_role")
        expect_identical(
            nodes[[slot_name]]$role_id,
            .artifactMetabSerializedSlots[[slot_name]]
        )
        expect_match(nodes[[slot_name]]$byte_digest, "^[a-f0-9]{64}$")
    }
    expect_false(any(vapply(nodes, \(node) {
        identical(node$node_type, "immutable_r_role") &&
            identical(node$slot_name, ".whole_object")
    }, logical(1))))
})

test_that("metabolomics codec rejects drift and future versions before hydration", {
    before <- .metabCodecDaResult(.metabCodecAssay("LCMS_Pos"))
    bundle <- dehydrateMetabolomicsS4Artifact(
        before,
        .metabCodecLogicalKey(stage_id = "da", state_role = "results")
    )
    future <- bundle
    future$metadata$codec$version <- 2L
    future$payloads[[1L]] <- NULL
    expect_error(
        hydrateMetabolomicsS4Artifact(future),
        class = "multischolar_unsupported_artifact_codec_version"
    )

    role_future <- bundle
    role_future$metadata$slot_values$fit.eb$codec$version <- 2L
    role_future$metadata$semantic_digest <- artifactSemanticDigest(
        artifactDiaBundleSemanticInput(role_future$metadata)
    )
    expect_error(
        hydrateMetabolomicsS4Artifact(role_future),
        class = "multischolar_unsupported_artifact_codec_version"
    )

    corrupt <- bundle
    role <- corrupt$metadata$slot_values$fit.eb
    payload <- as.data.frame(corrupt$payloads[[role$payload_key]])
    physical_name <- corrupt$metadata$payloads[[role$payload_key]]$
        columns[[2L]]$physical_name
    serialized_chunk <- jsonlite::base64_dec(payload[[physical_name]][[1L]])
    serialized_chunk[[length(serialized_chunk)]] <- as.raw(
        bitwXor(as.integer(serialized_chunk[[length(serialized_chunk)]]), 1L)
    )
    payload[[physical_name]][[1L]] <- jsonlite::base64_enc(serialized_chunk)
    corrupt$payloads[[role$payload_key]] <- payload
    expect_error(
        hydrateMetabolomicsS4Artifact(corrupt),
        class = "multischolar_artifact_hash_mismatch"
    )

    wrong_binding <- bundle
    wrong_binding$metadata$identity_binding$logical_key$generation_id <- "forged"
    wrong_binding$metadata$semantic_digest <- artifactSemanticDigest(
        artifactDiaBundleSemanticInput(wrong_binding$metadata)
    )
    expect_error(
        hydrateMetabolomicsS4Artifact(wrong_binding),
        class = "multischolar_metabolomics_codec_identity_mismatch"
    )
})

test_that("existing metabolomics RDS objects remain readable", {
    assay <- .metabCodecAssay(c("LCMS_Pos", "GCMS"))
    result <- .metabCodecDaResult(assay)
    assay_path <- tempfile(fileext = ".rds")
    result_path <- tempfile(fileext = ".rds")
    saveRDS(assay, assay_path, version = 3)
    saveRDS(result, result_path, version = 3)
    .metabCodecExpectExactSlots(assay, readRDS(assay_path))
    .metabCodecExpectDaSlots(result, readRDS(result_path))
})

.metabCodecWorkload <- function() {
    environment <- new.env(parent = .GlobalEnv)
    sys.source(
        .metabCodecRepoPath("tools", "profiling", "omics_workload_contract.R"),
        envir = environment
    )
    adapter_path <- .metabCodecRepoPath(
        "tools", "profiling", "omics_workload_metabolomics.R"
    )
    contract_path <- .metabCodecRepoPath(
        "tests", "testdata", "omics-parity", "metabolomics", "workloads",
        "mixed-public-ci-v1.json"
    )
    contract <- environment$omicsWorkloadReadContract(contract_path)
    adapter <- environment$omicsWorkloadLoadAdapter(adapter_path, contract)
    do.call(RNGkind, as.list(unlist(contract$rng$kind, use.names = FALSE)))
    set.seed(as.integer(contract$rng$seed))
    prepared <- adapter$prepare(list(
        contract = contract,
        run_dir = withr::local_tempdir(pattern = "metab-codec-workload-")
    ))
    payload_digest <- digest::digest(
        file = prepared$payload_path,
        algo = "sha256",
        serialize = FALSE
    )
    truth_digest <- digest::digest(
        file = prepared$truth_path,
        algo = "sha256",
        serialize = FALSE
    )
    .metabCodecAssertWorkloadDigests(payload_digest, truth_digest, contract)
    data <- utils::read.delim(
        prepared$payload_path,
        check.names = FALSE,
        stringsAsFactors = FALSE,
        na.strings = "NA"
    )
    truth <- jsonlite::read_json(prepared$truth_path, simplifyVector = TRUE)
    assays <- lapply(unlist(truth$assays, use.names = FALSE), \(assay_name) {
        assay <- data[data$assay == assay_name, , drop = FALSE]
        assay$assay <- NULL
        row.names(assay) <- NULL
        assay
    })
    names(assays) <- unlist(truth$assays, use.names = FALSE)
    design <- data.frame(
        sample_id = unlist(truth$sample_ids, use.names = FALSE),
        group = unlist(truth$group_assignments, use.names = FALSE),
        batch = unlist(truth$batch_assignments, use.names = FALSE),
        stringsAsFactors = FALSE
    )
    evidence <- list(
        workload_id = contract$workload_id,
        workload_digest = environment$omicsWorkloadDigest(contract),
        payload_sha256 = payload_digest,
        truth_sha256 = truth_digest,
        evidence_class = "generated_scaling_not_biological_validation"
    )
    list(
        object = createMetaboliteAssayData(
            metabolite_data = assays,
            design_matrix = design,
            metabolite_id_column = "feature_id",
            annotation_id_column = "annotation",
            sample_id = "sample_id",
            group_id = "group",
            technical_replicate_id = NA_character_,
            database_identifier_type = "Synthetic",
            internal_standard_regex = "^METAB_IS_",
            args = list(artifact_codec_evidence = evidence)
        ),
        evidence = evidence,
        contract = contract
    )
}

.metabCodecAssertWorkloadDigests <- function(
    payload_digest,
    truth_digest,
    contract
) {
    valid <- identical(
        payload_digest,
        contract$expected_digests$payload_sha256
    ) && identical(
        truth_digest,
        contract$expected_digests$truth_sha256
    )
    if (!isTRUE(valid)) {
        rlang::abort(
            "frozen metabolomics workload or truth digest drifted",
            class = "multischolar_metabolomics_workload_drift"
        )
    }
    invisible(TRUE)
}

test_that("frozen public workload round-trips with exact digest evidence", {
    workload <- .metabCodecWorkload()
    bundle <- dehydrateMetabolomicsS4Artifact(
        workload$object,
        .metabCodecLogicalKey(generation_id = "mixed-public-ci-v1")
    )
    restored <- hydrateMetabolomicsS4Artifact(
        .metabCodecWriteReadPayloads(bundle)
    )
    .metabCodecExpectExactSlots(workload$object, restored)
    expect_identical(
        restored@args$artifact_codec_evidence,
        workload$evidence
    )
    expect_identical(names(restored@metabolite_data), c(
        "LCMS_Pos", "LCMS_Neg", "GCMS"
    ))
    expect_identical(
        vapply(restored@metabolite_data, nrow, integer(1)),
        c(LCMS_Pos = 24L, LCMS_Neg = 24L, GCMS = 24L)
    )
    assay_nodes <- bundle$metadata$slot_values$metabolite_data$values
    assay_payload_keys <- vapply(
        assay_nodes,
        `[[`,
        character(1),
        "payload_key"
    )
    expect_identical(anyDuplicated(assay_payload_keys), 0L)
    expect_true(all(vapply(
        bundle$metadata$payloads[assay_payload_keys],
        \(metadata) metadata$writer_settings$chunk_size <= 65536L,
        logical(1)
    )))

    drifted <- workload$contract
    drifted$expected_digests$truth_sha256 <- strrep("0", 64L)
    expect_error(
        .metabCodecAssertWorkloadDigests(
            workload$evidence$payload_sha256,
            workload$evidence$truth_sha256,
            drifted
        ),
        class = "multischolar_metabolomics_workload_drift"
    )
})

test_that("metabolomics codecs are registered once without class changes", {
    catalogue <- artifactS4CodecCatalogue()
    ids <- names(artifactMetabolomicsCodecDeclarations())
    expect_true(all(ids %in% names(catalogue$codecs)))
    expect_identical(anyDuplicated(names(catalogue$codecs)), 0L)
    expect_identical(
        unname(vapply(catalogue$codecs[ids], `[[`, character(1), "class_name")),
        c("MetaboliteAssayData", "MetabolomicsDifferentialAbundanceResults")
    )
    description <- read.dcf(.metabCodecRepoPath("DESCRIPTION"))
    collate <- strsplit(description[1L, "Collate"], "[[:space:]]+")[[1L]]
    collate <- gsub("^'|'$", "", collate[nzchar(collate)])
    expect_identical(sum(collate == "utils_artifact_metab_codecs.R"), 1L)
    expect_lt(
        match("utils_artifact_s4_codecs.R", collate),
        match("utils_artifact_metab_codecs.R", collate)
    )
    source <- paste(readLines(
        .metabCodecRepoPath("R", "utils_artifact_metab_codecs.R"),
        warn = FALSE
    ), collapse = "\n")
    expect_false(grepl("saveRDS|readRDS", source))
    expect_false(grepl("serialize\\(value", source))
})
