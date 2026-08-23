.protNondiaCodecRepoPath <- function(...) {
    file.path(normalizePath(testthat::test_path("..", ".."), mustWork = TRUE), ...)
}

.protNondiaCodecManifest <- function() {
    jsonlite::read_json(
        .protNondiaCodecRepoPath(
            "tests", "testdata", "omics-parity", "proteomics", "manifest.json"
        ),
        simplifyVector = FALSE
    )
}

.protNondiaCodecImporter <- function(format) {
    switch(
        format,
        maxquant = importMaxQuantData,
        fragpipe = importFragPipeData,
        pd_tmt = importProteomeDiscovererTMTData
    )
}

.protNondiaCodecBuildObject <- function(scenario) {
    importer <- .protNondiaCodecImporter(scenario$format)
    imported <- suppressMessages(importer(
        .protNondiaCodecRepoPath(scenario$fixture_path)
    ))
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
    workflow_type <- if (identical(scenario$format, "pd_tmt")) "TMT" else "LFQ"
    manager <- new.env(parent = emptyenv())
    manager$saveState <- \(state_name, s4_data_object, ...) {
        manager$object <- s4_data_object
        invisible(state_name)
    }
    workflow <- list2env(
        list(
            data_cln = data,
            column_mapping = mapping,
            design_matrix = design,
            config_list = list(
                globalParameters = list(workflow_type = workflow_type)
            ),
            state_manager = manager
        ),
        parent = emptyenv()
    )
    suppressMessages(buildProtDesignStateCheckpoint(
        workflowData = workflow,
        workflowType = workflow_type,
        actionLabel = "non-DIA codec fixture",
        validateColumnMapping = TRUE
    ))
    manager$object
}

.protNondiaCodecNormalize <- function(value) {
    invisible(utils::capture.output(
        normalized <- suppressMessages(normaliseBetweenSamples(
            value,
            normalisation_method = "cyclicloess"
        ))
    ))
    normalized
}

.protNondiaCodecDecorate <- function(value, state_role) {
    sample_columns <- setdiff(
        names(value@protein_quant_table),
        value@protein_id_column
    )
    value@protein_quant_table[[sample_columns[[1L]]]][seq_len(4L)] <- c(
        NA_real_, NaN, Inf, -Inf
    )
    value@protein_id_table$review_status <- factor(
        rep(c("reviewed", "pending"), length.out = nrow(value@protein_id_table)),
        levels = c("pending", "reviewed")
    )
    value@design_matrix$group <- factor(
        value@design_matrix$group,
        levels = c("KO", "WT")
    )
    value@args$codec_fixture <- list(
        state_role = state_role,
        mapping = data.frame(
            source = c("protein", "sample"),
            target = c(value@protein_id_column, value@sample_id),
            stringsAsFactors = FALSE
        ),
        matrix = matrix(
            c(1, NA_real_, NaN, Inf, -Inf, 6),
            nrow = 2L,
            dimnames = list(c("first", "second"), c("a", "b", "c"))
        ),
        named_specials = c(
            missing = NA_real_,
            not_a_number = NaN,
            positive_infinity = Inf,
            negative_infinity = -Inf
        ),
        large_vector = seq_len(20000L)
    )
    expect_identical(methods::validObject(value, test = TRUE), TRUE)
    value
}

.protNondiaCodecWriteRead <- function(bundle) {
    payloads <- Map(\(payload, metadata, index) {
        encoded <- structure(
            list(payload = payload, metadata = metadata),
            class = c("MultiScholaRArtifactRectangular", "list")
        )
        path <- tempfile(
            sprintf("prot-nondia-codec-%03d-", index),
            fileext = ".parquet"
        )
        do.call(
            arrow::write_parquet,
            c(list(x = payload, sink = path), artifactParquetWriteArgs(encoded))
        )
        arrow::read_parquet(path, as_data_frame = FALSE)
    }, bundle$payloads, bundle$metadata$payloads, seq_along(bundle$payloads))
    names(payloads) <- names(bundle$payloads)
    bundle$payloads <- payloads
    bundle
}

.protNondiaCodecMetadataJsonRoundTrip <- function(bundle) {
    metadata_json <- jsonlite::toJSON(
        bundle$metadata,
        auto_unbox = TRUE,
        null = "null",
        na = "string",
        digits = NA
    )
    bundle$metadata <- jsonlite::fromJSON(
        metadata_json,
        simplifyVector = TRUE,
        simplifyDataFrame = FALSE,
        simplifyMatrix = FALSE
    )
    bundle
}

.expectProtNondiaCodecExact <- function(before, after) {
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

test_that("non-DIA codec roles cover only scientifically supported protein tuples", {
    roles <- artifactProteomicsNonDiaCodecRoles()
    expect_setequal(names(roles), c(
        "proteomics.maxquant.protein.lfq.v1",
        "proteomics.fragpipe.protein.lfq.v1",
        "proteomics.pd_tmt.protein.tmt.v1"
    ))
    expect_true(all(vapply(
        roles,
        \(role) identical(role$data_level, "protein") &&
            identical(role$s4_role, "protein_quantitative_state") &&
            identical(role$class_name, "ProteinQuantitativeData"),
        logical(1)
    )))
    expect_true(all(vapply(
        roles,
        \(role) identical(role$state_roles, .artifactProteomicsNonDiaStateRoles),
        logical(1)
    )))
    declarations <- artifactProteomicsNonDiaCodecDeclarations()
    catalogue <- artifactS4CodecCatalogue()
    expect_true(all(names(declarations) %in% names(catalogue$codecs)))
    expect_true(all(names(artifactDiaCodecDeclarations()) %in% names(catalogue$codecs)))
    expect_false(any(grepl("spectronaut", names(declarations), fixed = TRUE)))
    expect_false(any(vapply(
        declarations,
        \(codec) identical(codec$class_name, "PeptideQuantitativeData"),
        logical(1)
    )))
})

test_that("supported non-DIA protein states round-trip exactly by tuple and role", {
    testthat::skip_if_not_installed("arrow")
    manifest <- .protNondiaCodecManifest()
    for (scenario in manifest$fixture_scenarios) {
        initial <- .protNondiaCodecBuildObject(scenario)
        normalized <- .protNondiaCodecNormalize(initial)
        role <- artifactProteomicsNonDiaCodecRole(
            scenario$capability_id,
            "protein_s4_initial"
        )
        for (state_role in role$state_roles) {
            source <- if (identical(state_role, "protein_s4_initial")) {
                initial
            } else {
                normalized
            }
            before <- .protNondiaCodecDecorate(source, state_role)
            original <- unserialize(serialize(before, NULL))
            bundle <- dehydrateProteomicsNonDiaS4Artifact(
                before,
                scenario$capability_id,
                state_role
            )
            expect_identical(
                bundle$metadata$codec,
                list(id = role$codec_id, version = role$codec_version)
            )
            expect_identical(bundle$metadata$class_name, role$class_name)
            expect_gt(length(bundle$payloads), 3L)
            expect_false(any(vapply(bundle$payloads, is.raw, logical(1))))
            bundle <- .protNondiaCodecMetadataJsonRoundTrip(bundle)
            bundle <- .protNondiaCodecWriteRead(bundle)
            restored <- hydrateProteomicsNonDiaS4Artifact(
                bundle,
                scenario$capability_id,
                state_role
            )
            .expectProtNondiaCodecExact(before, restored)
            .expectProtNondiaCodecExact(original, before)
        }
    }
})

test_that("unsupported tuple, peptide, version, and state roles fail closed", {
    cases <- list(
        list(
            capability_id = "proteomics.spectronaut.protein.lfq.v1",
            state_role = "protein_s4_initial",
            version = "1.0.0"
        ),
        list(
            capability_id = "proteomics.spectronaut.peptide.lfq.v1",
            state_role = "raw_data_s4",
            version = "1.0.0"
        ),
        list(
            capability_id = "proteomics.maxquant.protein.lfq.v1",
            state_role = "raw_data_s4",
            version = "1.0.0"
        ),
        list(
            capability_id = "proteomics.maxquant.protein.lfq.v1",
            state_role = "protein_s4_initial",
            version = "2.0.0"
        ),
        list(
            capability_id = NA_character_,
            state_role = "protein_s4_initial",
            version = "1.0.0"
        )
    )
    for (case in cases) {
        error <- rlang::catch_cnd(artifactProteomicsNonDiaCodecRole(
            case$capability_id,
            case$state_role,
            case$version
        ))
        expect_s3_class(error, "multischolar_unsupported_proteomics_codec_role")
        expect_identical(error$capability_id, case$capability_id)
        expect_identical(error$state_role, case$state_role)
        expect_true(workflowCapabilityScalarString(error$remediation))
    }
})

test_that("unsupported non-DIA S4 shapes report lane, slot, and remediation", {
    scenario <- .protNondiaCodecManifest()$fixture_scenarios[[1L]]
    valid <- .protNondiaCodecBuildObject(scenario)

    reordered <- valid
    columns <- names(reordered@protein_quant_table)
    reordered@protein_quant_table <- reordered@protein_quant_table[
        c(columns[[1L]], rev(columns[-1L]))
    ]
    error <- rlang::catch_cnd(dehydrateProteomicsNonDiaS4Artifact(
        reordered,
        scenario$capability_id,
        "protein_s4_initial"
    ))
    expect_s3_class(error, "multischolar_proteomics_codec_shape_mismatch")
    expect_identical(error$capability_id, scenario$capability_id)
    expect_identical(error$state_role, "protein_s4_initial")
    expect_identical(error$slot_name, "protein_quant_table")
    expect_true(workflowCapabilityScalarString(error$remediation))

    missing_workflow <- valid
    missing_workflow@args$globalParameters$workflow_type <- NULL
    error <- rlang::catch_cnd(dehydrateProteomicsNonDiaS4Artifact(
        missing_workflow,
        scenario$capability_id,
        "protein_s4_initial"
    ))
    expect_identical(error$slot_name, "args$globalParameters$workflow_type")

    missing_provenance <- valid
    missing_provenance@protein_id_table <- missing_provenance@protein_id_table[-1L, ]
    error <- rlang::catch_cnd(dehydrateProteomicsNonDiaS4Artifact(
        missing_provenance,
        scenario$capability_id,
        "protein_s4_initial"
    ))
    expect_identical(
        error$slot_name,
        paste0("protein_id_table$", valid@protein_id_column)
    )

    opaque <- valid
    opaque@args$opaque_full_object_rds <- serialize(valid, NULL)
    error <- rlang::catch_cnd(dehydrateProteomicsNonDiaS4Artifact(
        opaque,
        scenario$capability_id,
        "protein_s4_initial"
    ))
    expect_s3_class(error, "multischolar_proteomics_codec_shape_mismatch")
    expect_match(error$slot_name, "opaque_full_object_rds", fixed = TRUE)

    peptide <- PeptideQuantitativeDataDiann(
        peptide_data = data.frame(
            Protein.Group = c("PG1", "PG1"),
            Stripped.Sequence = c("PEPTIDE", "PEPTIDE"),
            Run = c("S1", "S2"),
            Q.Value = c(0.01, 0.02),
            Precursor.Quantity = c(10, 20),
            Precursor.Normalised = c(10, 20),
            stringsAsFactors = FALSE
        ),
        design_matrix = data.frame(
            Run = c("S1", "S2"),
            group = c("WT", "KO"),
            replicates = c("R1", "R1"),
            stringsAsFactors = FALSE
        ),
        args = list(globalParameters = list(workflow_type = "LFQ"))
    )
    error <- rlang::catch_cnd(dehydrateProteomicsNonDiaS4Artifact(
        peptide,
        scenario$capability_id,
        "protein_s4_initial"
    ))
    expect_identical(error$slot_name, ".class")
    expect_match(error$remediation, "protein-level", fixed = TRUE)
})

test_that("codec metadata cannot be borrowed across supported tuples", {
    manifest <- .protNondiaCodecManifest()
    maxquant <- manifest$fixture_scenarios[[1L]]
    fragpipe <- manifest$fixture_scenarios[[2L]]
    object <- .protNondiaCodecBuildObject(maxquant)
    bundle <- dehydrateProteomicsNonDiaS4Artifact(
        object,
        maxquant$capability_id,
        "protein_s4_initial"
    )
    expect_error(
        hydrateProteomicsNonDiaS4Artifact(
            bundle,
            fragpipe$capability_id,
            "protein_s4_initial"
        ),
        class = "multischolar_unsupported_artifact_codec_version"
    )
})

test_that("legacy RDS readback and memory construction remain unchanged", {
    for (scenario in .protNondiaCodecManifest()$fixture_scenarios) {
        before <- .protNondiaCodecBuildObject(scenario)
        path <- tempfile("prot-nondia-legacy-", fileext = ".rds")
        saveRDS(before, path)
        restored <- readRDS(path)
        .expectProtNondiaCodecExact(before, restored)
        expect_identical(
            methods::validObject(.protNondiaCodecBuildObject(scenario), test = TRUE),
            TRUE
        )
    }
})
