# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

.artifactLipidomicsCapabilityVersion <- "1.0.0"

artifactLipidomicsCodecDeclarations <- function() {
    list(
        "multischolar.s4.lipidomics_assay_data" = list(
            codec_id = "multischolar.s4.lipidomics_assay_data",
            codec_version = 1L,
            class_name = "LipidomicsAssayData",
            payload_schema_id = "multischolar.rectangular",
            payload_schema_version = 1L
        ),
        "multischolar.s4.lipidomics_da_results" = list(
            codec_id = "multischolar.s4.lipidomics_da_results",
            codec_version = 1L,
            class_name = "LipidomicsDifferentialAbundanceResults",
            payload_schema_id = "multischolar.rectangular",
            payload_schema_version = 1L
        )
    )
}

artifactLipidomicsImportQuery <- function() {
    operation_id <- "lipidomics.lipidsearch.import.assay_manifest.v1"
    stats::setNames(list(list(
        operation_id = operation_id,
        state_role = "assay_manifest",
        projections = c(
            "assay_name", "assay_order", "artifact_role", "feature_count",
            "sample_count", "lipid_id_column", "source_digest"
        ),
        filters = list(assay_name = list(
            column = "assay_name",
            type = "character",
            operators = c("equal", "in")
        )),
        order_by = "assay_order",
        max_rows = 1000L,
        max_bytes = 16L * 1024L * 1024L
    )), operation_id)
}

artifactLipidomicsPerformanceThresholds <- function() {
    c(
        maximum_peak_tree_rss_bytes = 1132247040,
        maximum_retained_tree_rss_bytes = 1132247040,
        maximum_elapsed_seconds = 18.17375,
        maximum_final_disk_bytes = 3174,
        maximum_final_file_count = 5,
        maximum_bounded_query_p95_seconds = 0.0012500000000018052
    )
}

artifactLipidomicsWorkflowDescriptor <- function() {
    owner <- "lipidomics.lipidsearch.lipid.standard.v1"
    codec_id <- "multischolar.s4.lipidomics_assay_data"
    codecs <- artifactLipidomicsCodecDeclarations()[codec_id]
    queries <- artifactLipidomicsImportQuery()
    query_id <- names(queries)[[1L]]
    newArtifactWorkflowDescriptor(
        descriptor_id = owner,
        descriptor_version = .artifactLipidomicsCapabilityVersion,
        identity = workflowCapabilityIdentity(
            "lipidomics",
            "lipidomics.gui",
            "lipidomics_standard",
            "lipid_standard",
            "lipidsearch",
            "lipid",
            "not_recorded"
        ),
        stages = list(
            import = list(
                stage_id = "import",
                state_roles = c(
                    "assay_payloads", "assay_manifest", "column_mapping",
                    "source_manifest"
                ),
                codec_ids = codec_id,
                query_operation_ids = query_id,
                maximum_rollout = "evict"
            ),
            design = list(
                stage_id = "design",
                state_roles = c(
                    "cleaned_assay_payloads", "design_matrix", "contrasts",
                    "args", "column_mapping", "assay_alignment",
                    "raw_s4_dependencies", "lipid_raw_data_s4"
                ),
                codec_ids = codec_id,
                query_operation_ids = character(),
                maximum_rollout = "evict"
            )
        ),
        codecs = codecs,
        queries = queries,
        fixtures = list(
            owner_id = owner,
            fixture_ids = c(
                "tests/testdata/e2e/lipid_canonical/lipidsearch_lcms_pos.txt",
                "tests/testdata/e2e/lipid_canonical/lipidsearch_lcms_neg.txt",
                "tests/testdata/e2e/lipid_canonical/lipidsearch_gcms.txt"
            )
        ),
        scientific_oracle = list(
            owner_id = owner,
            oracle_id = "tests/testdata/omics-parity/lipidomics/memory-oracle.json",
            oracle_version = "1.0.0",
            tolerances = c(
                import_absolute = 1e-8,
                import_relative = 1e-8,
                differential_abundance_absolute = 1e-7
            )
        ),
        compatibility_products = list(
            owner_id = owner,
            product_ids = c(
                "data_cln_*.tab", "assay_manifest.txt", "column_mapping.json",
                "lipidomics_import_summary.tsv", "design_matrix.tab",
                "contrast_strings.tab", "config.ini", "manifest.json"
            )
        ),
        evidence = list(
            owner_id = owner,
            inventory_ids = c(
                "tests/testdata/omics-capabilities.json",
                "tests/testdata/omics-parity/lipidomics/manifest.json",
                "tests/testdata/omics-parity/lipidomics/resources/ls-pos-memory-v1.json",
                "tests/testdata/omics-parity/lipidomics/resources/ls-neg-memory-v1.json",
                "tests/testdata/omics-parity/lipidomics/resources/ls-gc-memory-v1.json",
                paste0(
                    "tests/testdata/omics-parity/lipidomics/workloads/",
                    "mixed-public-ci-v1.json"
                )
            ),
            codec_ids = codec_id,
            stage_ids = c("import", "design"),
            lifecycle_ids = c("OMICS-ART-038", "OMICS-ART-058"),
            performance_thresholds = artifactLipidomicsPerformanceThresholds()
        ),
        migration = list(
            owner_id = owner,
            strategy_id = "explicit_lipidsearch_lipidomics_dual_write",
            from_backend = "memory",
            to_backend = "artifact"
        ),
        rollback = list(
            owner_id = owner,
            strategy_id = "force_memory_ignore_lipidomics_generations",
            target_backend = "memory"
        ),
        certification = list(
            owner_id = owner,
            status = "evict",
            auto_eligible = FALSE
        )
    )
}
