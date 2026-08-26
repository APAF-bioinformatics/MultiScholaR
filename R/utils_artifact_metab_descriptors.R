# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

.artifactMetabolomicsCapabilityVersion <- "1.0.0"

artifactMetabolomicsCodecDeclarations <- function() {
    list(
        "multischolar.s4.metabolite_assay_data" = list(
            codec_id = "multischolar.s4.metabolite_assay_data",
            codec_version = 1L,
            class_name = "MetaboliteAssayData",
            payload_schema_id = "multischolar.rectangular",
            payload_schema_version = 1L
        ),
        "multischolar.s4.metabolomics_da_results" = list(
            codec_id = "multischolar.s4.metabolomics_da_results",
            codec_version = 1L,
            class_name = "MetabolomicsDifferentialAbundanceResults",
            payload_schema_id = "multischolar.rectangular",
            payload_schema_version = 1L
        )
    )
}

artifactMetabolomicsImportQuery <- function() {
    operation_id <- "metabolomics.custom.import.assay_manifest.v1"
    stats::setNames(list(list(
        operation_id = operation_id,
        state_role = "assay_manifest",
        projections = c(
            "assay_name", "assay_order", "artifact_role", "feature_count",
            "sample_count", "feature_id_column", "source_digest"
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

artifactMetabolomicsPerformanceThresholds <- function() {
    c(
        maximum_peak_tree_rss_bytes = 1149864960,
        maximum_retained_tree_rss_bytes = 1134709760,
        maximum_elapsed_seconds = 8.573,
        maximum_final_disk_bytes = 3133,
        maximum_final_file_count = 5,
        maximum_bounded_query_p95_seconds = 0.447
    )
}

artifactMetabolomicsWorkflowDescriptor <- function() {
    owner <- "metabolomics.custom.metabolite.standard.v1"
    codec_id <- "multischolar.s4.metabolite_assay_data"
    codecs <- artifactMetabolomicsCodecDeclarations()[codec_id]
    queries <- artifactMetabolomicsImportQuery()
    query_id <- names(queries)[[1L]]
    newArtifactWorkflowDescriptor(
        descriptor_id = owner,
        descriptor_version = .artifactMetabolomicsCapabilityVersion,
        identity = workflowCapabilityIdentity(
            "metabolomics",
            "metabolomics.gui",
            "metabolomics_standard",
            "metab_standard",
            "custom",
            "metabolite",
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
                    "raw_s4_dependencies", "metab_raw_data_s4"
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
                "tests/testdata/e2e/metab_lc/lcms_pos_features.tsv",
                "tests/testdata/e2e/metab_lc/lcms_neg_features.tsv",
                "tests/testdata/e2e/metab_gc/gcms_features.tsv",
                "tests/testdata/e2e/metab_combined/combined_lcms_features.tsv",
                "tests/testdata/e2e/metab_combined/combined_gcms_features.tsv"
            )
        ),
        scientific_oracle = list(
            owner_id = owner,
            oracle_id = "tests/testdata/omics-parity/metabolomics/memory-oracle.json",
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
                "metabolomics_import_summary.tsv", "design_matrix.tab",
                "contrast_strings.tab", "config.ini", "manifest.json"
            )
        ),
        evidence = list(
            owner_id = owner,
            inventory_ids = c(
                "tests/testdata/omics-capabilities.json",
                "tests/testdata/omics-parity/metabolomics/manifest.json",
                paste0(
                    "tests/testdata/omics-parity/metabolomics/resources/",
                    "mixed-memory-baseline-v1.json"
                ),
                paste0(
                    "tests/testdata/omics-parity/metabolomics/workloads/",
                    "mixed-public-ci-v1.json"
                )
            ),
            codec_ids = codec_id,
            stage_ids = c("import", "design"),
            lifecycle_ids = c("OMICS-ART-029", "OMICS-ART-058"),
            performance_thresholds = artifactMetabolomicsPerformanceThresholds()
        ),
        migration = list(
            owner_id = owner,
            strategy_id = "explicit_custom_metabolomics_dual_write",
            from_backend = "memory",
            to_backend = "artifact"
        ),
        rollback = list(
            owner_id = owner,
            strategy_id = "force_memory_ignore_metabolomics_generations",
            target_backend = "memory"
        ),
        certification = list(
            owner_id = owner,
            status = "evict",
            auto_eligible = FALSE
        )
    )
}
