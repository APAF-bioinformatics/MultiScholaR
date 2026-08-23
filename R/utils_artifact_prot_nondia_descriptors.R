# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

#' Define supported non-DIA proteomics descriptor specifications
#'
#' @return A named list of immutable workflow, evidence, and resource contracts.
#' @noRd
artifactProteomicsNonDiaDescriptorSpecs <- function() {
    list(
        "proteomics.maxquant.protein.lfq.v1" = list(
            workflow_type = "LFQ",
            workflow_slug = "prot_lfq",
            input_format = "maxquant",
            codec_id = paste0(
                "multischolar.s4.protein_quantitative_data.",
                "maxquant.protein.lfq"
            ),
            fixture_ids = c(
                "tests/testdata/e2e/prot_lfq/proteinGroups.txt",
                "tests/testdata/e2e/prot_lfq/design.tsv"
            ),
            resource_id = paste0(
                "tests/testdata/omics-parity/proteomics/resources/",
                "maxquant-memory-baseline-v1.json"
            ),
            workload_id = paste0(
                "tests/testdata/omics-parity/proteomics/workloads/",
                "maxquant-protein-lfq-public-v1.json"
            ),
            projections = c("Protein.Ids", "Gene.names", "Run", "Intensity"),
            thresholds = c(
                maximum_peak_tree_rss_bytes = 1167503360,
                maximum_retained_tree_rss_bytes = 1167503360,
                maximum_elapsed_seconds = 18.698,
                maximum_peak_disk_bytes = 3400,
                maximum_final_disk_bytes = 3400,
                maximum_final_file_count = 5,
                maximum_bounded_query_p95_seconds = 0.093
            )
        ),
        "proteomics.fragpipe.protein.lfq.v1" = list(
            workflow_type = "LFQ",
            workflow_slug = "prot_lfq_fragpipe",
            input_format = "fragpipe",
            codec_id = paste0(
                "multischolar.s4.protein_quantitative_data.",
                "fragpipe.protein.lfq"
            ),
            fixture_ids = c(
                "tests/testdata/e2e/prot_lfq/seed_combined_protein.tsv",
                "tests/testdata/e2e/prot_lfq/design.tsv"
            ),
            resource_id = paste0(
                "tests/testdata/omics-parity/proteomics/resources/",
                "fragpipe-memory-baseline-v1.json"
            ),
            workload_id = paste0(
                "tests/testdata/omics-parity/proteomics/workloads/",
                "fragpipe-protein-lfq-public-v1.json"
            ),
            projections = c("Protein.Ids", "Run", "Intensity"),
            thresholds = c(
                maximum_peak_tree_rss_bytes = 1164764160,
                maximum_retained_tree_rss_bytes = 1164764160,
                maximum_elapsed_seconds = 18.697,
                maximum_peak_disk_bytes = 4139,
                maximum_final_disk_bytes = 4139,
                maximum_final_file_count = 5,
                maximum_bounded_query_p95_seconds = 0.103
            )
        ),
        "proteomics.pd_tmt.protein.tmt.v1" = list(
            workflow_type = "TMT",
            workflow_slug = "prot_tmt",
            input_format = "pd_tmt",
            codec_id = paste0(
                "multischolar.s4.protein_quantitative_data.",
                "pd_tmt.protein.tmt"
            ),
            fixture_ids = c(
                "tests/testdata/e2e/prot_tmt/pd_tmt_peptides.tsv",
                "tests/testdata/e2e/prot_tmt/design.tsv"
            ),
            resource_id = paste0(
                "tests/testdata/omics-parity/proteomics/resources/",
                "pd_tmt-memory-baseline-v1.json"
            ),
            workload_id = paste0(
                "tests/testdata/omics-parity/proteomics/workloads/",
                "pd-tmt-protein-public-v1.json"
            ),
            projections = c(
                "Annotated.Sequence", "Protein.Ids", "Run", "Abundance"
            ),
            thresholds = c(
                maximum_peak_tree_rss_bytes = 1169105920,
                maximum_retained_tree_rss_bytes = 1169105920,
                maximum_elapsed_seconds = 19.225,
                maximum_peak_disk_bytes = 4134,
                maximum_final_disk_bytes = 4134,
                maximum_final_file_count = 5,
                maximum_bounded_query_p95_seconds = 0.147
            )
        )
    )
}

#' Build the bounded import query contract for a non-DIA tuple
#'
#' @param owner Descriptor identifier that owns the query contract.
#' @param spec Supported tuple specification.
#'
#' @return A named list containing one bounded query declaration.
#' @noRd
artifactProteomicsNonDiaImportQuery <- function(owner, spec) {
    operation_id <- paste0(owner, ".import.preview.v1")
    projections <- spec$projections
    filters <- list(
        protein = list(
            column = "Protein.Ids",
            type = "character",
            operators = c("equal", "in")
        ),
        run = list(
            column = "Run",
            type = "character",
            operators = c("equal", "in")
        )
    )
    stats::setNames(list(list(
        operation_id = operation_id,
        state_role = "canonical_data",
        projections = projections,
        filters = filters,
        order_by = c("Protein.Ids", "Run"),
        max_rows = 10000L,
        max_bytes = 64L * 1024L * 1024L
    )), operation_id)
}

#' Construct an immutable non-DIA proteomics workflow descriptor
#'
#' @param capability_id Exact supported workflow capability identifier.
#'
#' @return A validated `ArtifactWorkflowDescriptor` object.
#' @noRd
artifactProteomicsNonDiaWorkflowDescriptor <- function(capability_id) {
    specs <- artifactProteomicsNonDiaDescriptorSpecs()
    spec <- specs[[capability_id]]
    if (is.null(spec)) {
        artifactDescriptorAbort(
            "non-DIA proteomics descriptor has no supported exact tuple",
            "multischolar_missing_artifact_descriptor",
            capability_id = capability_id
        )
    }
    codecs <- artifactProteomicsNonDiaCodecDeclarations()
    codec <- codecs[spec$codec_id]
    queries <- artifactProteomicsNonDiaImportQuery(capability_id, spec)
    query_id <- names(queries)[[1L]]
    newArtifactWorkflowDescriptor(
        descriptor_id = capability_id,
        descriptor_version = .artifactProteomicsNonDiaCapabilityVersion,
        identity = workflowCapabilityIdentity(
            "proteomics",
            "proteomics.gui",
            spec$workflow_type,
            spec$workflow_slug,
            spec$input_format,
            "protein",
            "not_recorded"
        ),
        stages = list(
            import = list(
                stage_id = "import",
                state_roles = "canonical_data",
                codec_ids = spec$codec_id,
                query_operation_ids = query_id,
                maximum_rollout = "dual_write"
            ),
            design = list(
                stage_id = "design",
                state_roles = c(
                    "cleaned_data", "design_matrix", "contrasts", "args",
                    "annotations", "sequences", "protein_s4_initial"
                ),
                codec_ids = spec$codec_id,
                query_operation_ids = character(),
                maximum_rollout = "dual_write"
            )
        ),
        codecs = codec,
        queries = queries,
        fixtures = list(
            owner_id = capability_id,
            fixture_ids = spec$fixture_ids
        ),
        scientific_oracle = list(
            owner_id = capability_id,
            oracle_id = "tests/testdata/omics-parity/proteomics/memory-oracle.json",
            oracle_version = "1.0.0",
            tolerances = c(
                import_absolute = 1e-10,
                import_relative = 1e-12,
                normalization_absolute = 1e-8,
                differential_abundance_absolute = 1e-7
            )
        ),
        compatibility_products = list(
            owner_id = capability_id,
            product_ids = c(
                "data_cln.tab", "design_matrix.tab", "contrast_strings.tab",
                "config.ini", "aa_seq_tbl_final.RDS", "uniprot_dat_cln.RDS"
            )
        ),
        evidence = list(
            owner_id = capability_id,
            inventory_ids = c(
                "tests/testdata/omics-capabilities.json",
                "tests/testdata/omics-parity/proteomics/manifest.json",
                spec$resource_id,
                spec$workload_id
            ),
            codec_ids = spec$codec_id,
            stage_ids = c("import", "design"),
            lifecycle_ids = c("OMICS-ART-021", "OMICS-ART-058"),
            performance_thresholds = spec$thresholds
        ),
        migration = list(
            owner_id = capability_id,
            strategy_id = "explicit_tuple_dual_write",
            from_backend = "memory",
            to_backend = "artifact"
        ),
        rollback = list(
            owner_id = capability_id,
            strategy_id = "force_memory_ignore_tuple_generations",
            target_backend = "memory"
        ),
        certification = list(
            owner_id = capability_id,
            status = "dual_write",
            auto_eligible = FALSE
        )
    )
}

#' Construct all certified non-DIA proteomics workflow descriptors
#'
#' @return A named list of validated workflow descriptors.
#' @noRd
artifactProteomicsNonDiaWorkflowDescriptors <- function() {
    descriptors <- lapply(
        names(artifactProteomicsNonDiaDescriptorSpecs()),
        artifactProteomicsNonDiaWorkflowDescriptor
    )
    stats::setNames(
        descriptors,
        vapply(descriptors, `[[`, character(1), "descriptor_id")
    )
}
