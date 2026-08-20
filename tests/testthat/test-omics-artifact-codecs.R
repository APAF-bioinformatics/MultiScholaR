writeReadArtifactPayloads <- function(bundle) {
    payloads <- Map(function(payload, metadata, index) {
        encoded <- structure(
            list(payload = payload, metadata = metadata),
            class = c("MultiScholaRArtifactRectangular", "list")
        )
        path <- tempfile(sprintf("artifact-payload-%03d-", index), fileext = ".parquet")
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

expectExactS4Slots <- function(before, after) {
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

makeDiaCodecObjects <- function(branch = c("iq", "limpa")) {
    branch <- match.arg(branch)
    peptide_data <- data.frame(
        Protein.Group = rep(c("PG1", "PG2"), each = 2L),
        Protein.Ids = rep(c("P1;P1A", "P2"), each = 2L),
        Stripped.Sequence = rep(c("PEPTIDE_A", "PEPTIDE_B"), each = 2L),
        Run = rep(c("S1", "S2"), 2L),
        Q.Value = c(0.001, 0.002, 0.003, 0.004),
        Global.Q.Value = c(0.005, 0.006, 0.007, 0.008),
        Proteotypic = c(1L, 1L, 0L, 0L),
        Precursor.Quantity = c(100, NA_real_, 300, 400),
        Precursor.Normalised = c(10, 20, 30, 40),
        stringsAsFactors = FALSE
    )
    design <- data.frame(
        Run = c("S1", "S2"),
        group = factor(c("control", "case"), levels = c("case", "control")),
        replicates = c("R1", "R1"),
        stringsAsFactors = FALSE
    )
    peptide <- PeptideQuantitativeDataDiann(
        peptide_data,
        design,
        sample_id = "Run",
        group_id = "group",
        technical_replicate_id = "replicates",
        args = list(branch = branch)
    )
    peptide@peptide_matrix <- matrix(
        c(10, 30, 20, 40),
        nrow = 2L,
        dimnames = list(c("PG1%PEPTIDE_A", "PG2%PEPTIDE_B"), c("S1", "S2"))
    )
    peptide <- .ensurePeptideFeatureKeyMap(peptide, record_migration = FALSE)
    peptide@args$peptide_qc_audit <- list(
        active_digest = paste0(branch, "-audit-digest"),
        records = list(list(stage = "imputation", retained = 4L))
    )
    if (identical(branch, "iq")) {
        peptide@args$iq <- list(
            method = "maxLFQ",
            fitted = matrix(c(10, 30, 20, 40), nrow = 2L),
            summary = data.frame(protein = c("PG1", "PG2"), peptides = c(1L, 1L))
        )
    } else {
        peptide@args$limpa_dpc_results <- list(
            dpc_parameters = c(beta0 = -1.2, beta1 = 0.8),
            dpc_object = structure(
                list(dpc = c(-1.2, 0.8), fitted = matrix(c(0.1, 0.9), nrow = 1L)),
                class = "DPCFixture"
            ),
            slope_interpretation = "intensity-dependent"
        )
    }

    protein_table <- data.frame(
        Protein.Group = c("PG1", "PG2"),
        S1 = c(10, 30),
        S2 = c(20, 40),
        check.names = FALSE
    )
    if (identical(branch, "iq")) {
        protein_args <- list(
            branch = branch,
            peptide_feature_key_map = peptide@args$peptide_feature_key_map,
            peptide_qc_audit = peptide@args$peptide_qc_audit
        )
        protein_args$iq_results <- list(
            method = "maxLFQ",
            annotation = data.frame(Protein.Group = c("PG1", "PG2"), n = c(1L, 1L))
        )
        protein <- ProteinQuantitativeData(
            protein_quant_table = protein_table,
            protein_id_column = "Protein.Group",
            design_matrix = design,
            protein_id_table = data.frame(
                Protein.Group = c("PG1", "PG2"),
                Protein.Ids = c("P1;P1A", "P2"),
                stringsAsFactors = FALSE
            ),
            sample_id = "Run",
            group_id = "group",
            technical_replicate_id = "replicates",
            args = protein_args
        )
    } else {
        protein <- buildProtTestModeLimpaProteinObject(peptide, dpcSlope = 0.8)
    }
    list(peptide = peptide, protein = protein)
}

test_that("artifact refs are scoped, versioned, contained, and integrity checked", {
    project_root <- withr::local_tempdir()
    dir.create(file.path(project_root, "data", "proteomics"), recursive = TRUE)
    relative_path <- "data/proteomics/payload.parquet"
    payload_path <- file.path(project_root, relative_path)
    writeBin(charToRaw("immutable-payload"), payload_path)
    shape <- list(
        kind = "data.frame",
        rows = 4L,
        columns = 3L,
        payloads = 1L,
        bytes = unname(file.info(payload_path)$size)
    )
    semantic_digest <- artifactSemanticDigest(list(stage = "qc", rows = 4L))
    ref <- newArtifactRef(
        logical_key = list(
            project_id = "project-001",
            omic_type = "proteomics",
            workflow_slug = "diann",
            stage_id = "peptide_qc",
            state_role = "filtered_peptides",
            generation_id = "generation-001"
        ),
        relative_path = relative_path,
        codec_id = "multischolar.s4.peptide_quantitative_data.diann",
        payload_schema_id = "multischolar.parquet_table",
        shape = shape,
        semantic_digest = semantic_digest,
        byte_digest = artifactByteDigest(payload_path),
        status = "committed"
    )

    expect_s3_class(validateArtifactRef(ref), "MultiScholaRArtifactRef")
    expect_true(grepl("^art_[a-f0-9]{64}$", ref$artifact_id))
    expect_identical(ref$relative_path, relative_path)
    expect_identical(ref$hash_policy$semantic$digest, semantic_digest)
    expect_false(identical(
        ref$hash_policy$byte$digest,
        ref$hash_policy$semantic$digest
    ))
    expect_true(endsWith(ref$created_at, "Z"))
    expect_invisible(validateArtifactRefPayload(ref, project_root, shape))

    bad_version <- ref
    bad_version$schema_version <- 2L
    expect_error(
        validateArtifactRef(bad_version),
        class = "multischolar_unsupported_artifact_ref_version"
    )
    bad_shape <- shape
    bad_shape$rows <- 5L
    expect_error(
        validateArtifactRefPayload(ref, project_root, bad_shape),
        class = "multischolar_artifact_shape_mismatch"
    )
    writeBin(charToRaw("corrupt"), payload_path)
    expect_error(
        validateArtifactRefPayload(ref, project_root, shape),
        class = "multischolar_artifact_hash_mismatch"
    )
    expect_error(
        newArtifactRef(
            ref$logical_key,
            "../escape.parquet",
            ref$codec$id,
            ref$payload_schema$id,
            shape,
            semantic_digest
        ),
        class = "multischolar_invalid_relative_artifact_path"
    )
})

test_that("table codecs preserve exact R semantics through Parquet", {
    skip_if_not_installed("bit64")
    timestamp <- as.POSIXct(
        c("2020-01-01 00:00:00", NA, "2020-06-01 12:30:00", "2021-01-01 00:00:00"),
        tz = "Australia/Sydney"
    )
    value <- data.frame(
        numeric = c(1, NA_real_, NaN, Inf),
        negative_infinity = c(-Inf, 2, 3, 4),
        integer = c(1L, NA_integer_, 3L, 4L),
        logical = c(TRUE, FALSE, NA, TRUE),
        character = c("NA", NA_character_, "", "value"),
        factor = factor(c("b", NA, "a", "b"), levels = c("b", "a", "unused")),
        ordered = ordered(c("low", "high", NA, "low"), levels = c("low", "high")),
        date = as.Date(c("2020-01-01", NA, "2020-06-01", "2021-01-01")),
        timestamp = timestamp,
        integer64 = bit64::as.integer64(c("9007199254740993", NA, "3", "4")),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    names(value)[c(1L, 2L)] <- c(
        ".multischolar_row_order",
        ".multischolar_value_00000001"
    )
    names(value)[c(3L, 4L)] <- "duplicate"
    row.names(value) <- c("row-4", "row-3", "row-2", "row-1")
    attr(value[[1L]], "units") <- "intensity"

    encoded <- encodeArtifactTable(value)
    path <- tempfile(fileext = ".parquet")
    do.call(
        arrow::write_parquet,
        c(list(x = encoded$payload, sink = path), artifactParquetWriteArgs(encoded))
    )
    restored <- decodeArtifactRectangular(
        arrow::read_parquet(path, as_data_frame = FALSE),
        encoded$metadata
    )

    expect_identical(restored, value)
    expect_identical(encoded$metadata$logical_names, names(value))
    expect_false(anyDuplicated(vapply(
        encoded$metadata$physical_schema,
        `[[`,
        character(1),
        "name"
    )) > 0L)
    expect_true(all(c(
        "physical_schema", "stable_key", "schema_evolution", "writer_settings",
        "semantic_digest"
    ) %in% names(encoded$metadata)))
    expect_identical(encoded$metadata$schema_evolution$policy, "reject_unknown")
    expect_identical(encoded$metadata$writer_settings$append_policy, "immutable_new_generation")
    expect_identical(encoded$metadata$writer_settings$partitioning, "none")
})

test_that("table codecs preserve zero-row automatic row names", {
    value <- data.frame(
        artifact_role = character(),
        availability = character(),
        stringsAsFactors = FALSE
    )
    encoded <- encodeArtifactTable(value)
    restored <- decodeArtifactRectangular(encoded$payload, encoded$metadata)

    expect_identical(encoded$metadata$row_names$kind, "automatic")
    expect_identical(restored, value)
})

test_that("matrix codecs preserve dimensions, dimnames, storage, and non-finite values", {
    value <- matrix(
        c(1, NA_real_, NaN, Inf, -Inf, 6),
        nrow = 2L,
        dimnames = list(c("row-b", "row-a"), c("same", "same", "reserved"))
    )
    encoded <- encodeArtifactMatrix(value)
    restored <- decodeArtifactRectangular(encoded$payload, encoded$metadata)
    expect_identical(restored, value)
    expect_identical(encoded$metadata$matrix$storage_mode, "double")
    expect_identical(encoded$metadata$stable_key$kind, "artifact_row_order")
})

test_that("writer policy is explicit for representative table shapes", {
    fixtures <- list(
        wide = as.data.frame(matrix(seq_len(32L * 128L), nrow = 32L)),
        tall = data.frame(id = seq_len(5000L), group = rep(letters[1:4], length.out = 5000L)),
        sparse = data.frame(
            id = seq_len(2048L),
            value = replace(rep(NA_real_, 2048L), seq(1L, 2048L, by = 128L), 1)
        ),
        high_cardinality = data.frame(
            id = seq_len(2000L),
            label = sprintf("feature-%08d", seq_len(2000L))
        )
    )
    encoded <- lapply(names(fixtures), function(name) {
        encodeArtifactTable(fixtures[[name]], owner = name)
    })
    names(encoded) <- names(fixtures)

    for (name in names(encoded)) {
        item <- encoded[[name]]
        path <- tempfile(fileext = ".parquet")
        do.call(
            arrow::write_parquet,
            c(list(x = item$payload, sink = path), artifactParquetWriteArgs(item))
        )
        expect_identical(
            decodeArtifactRectangular(
                arrow::read_parquet(path, as_data_frame = FALSE),
                item$metadata
            ),
            fixtures[[name]],
            info = name
        )
        settings <- item$metadata$writer_settings
        expect_lte(settings$chunk_size, 65536L)
        expect_true(settings$write_statistics)
        expect_true(settings$compression %in% c("zstd", "snappy"))
        expect_identical(settings$partitioning, "none")
        expect_identical(settings$small_file_policy, "one_payload_file_per_generation")
    }
    expect_false(encoded$high_cardinality$metadata$writer_settings$use_dictionary)
    expect_identical(
        encoded$high_cardinality$metadata$writer_settings$dictionary_columns,
        character()
    )
})

test_that("DIA IQ and limpa peptide and protein objects hydrate exactly", {
    for (branch in c("iq", "limpa")) {
        objects <- makeDiaCodecObjects(branch)
        for (object_name in names(objects)) {
            before <- objects[[object_name]]
            bundle <- dehydrateDiaS4Artifact(before)
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
            restored <- hydrateDiaS4Artifact(writeReadArtifactPayloads(bundle))
            expectExactS4Slots(before, restored)
            expect_identical(restored@args, before@args)
            expect_identical(
                restored@args$peptide_qc_audit$active_digest,
                paste0(branch, "-audit-digest")
            )
            expect_identical(
                restored@args$peptide_feature_key_map,
                before@args$peptide_feature_key_map
            )
        }
    }
})

test_that("nested args externalize rectangular values and fail safely when unsupported", {
    objects <- makeDiaCodecObjects("limpa")
    objects$peptide@args$large_vector <- seq_len(20000L)
    bundle <- dehydrateDiaS4Artifact(objects$peptide, inline_limit_bytes = 65536L)
    restored <- hydrateDiaS4Artifact(bundle)
    expect_identical(restored@args$large_vector, seq_len(20000L))
    expect_gt(length(bundle$payloads), 2L)

    objects$peptide@args$unsafe_callback <- function(value) value
    expect_error(
        dehydrateDiaS4Artifact(objects$peptide),
        regexp = "PeptideQuantitativeData@args.*unsafe_callback",
        class = "multischolar_artifact_externalization_required"
    )
    oversized <- makeDiaCodecObjects("iq")$peptide
    oversized@args <- list(many_scalars = as.list(seq_len(500L)))
    oversized_bundle <- dehydrateDiaS4Artifact(
        oversized,
        inline_limit_bytes = 512L
    )
    expect_gt(length(oversized_bundle$payloads), 2L)
    expectExactS4Slots(oversized, hydrateDiaS4Artifact(oversized_bundle))

    unsafe_nested <- oversized
    unsafe_nested@args$many_scalars[[500L]] <- "/tmp/unsafe-artifact-path"
    expect_error(
        dehydrateDiaS4Artifact(unsafe_nested, inline_limit_bytes = 512L),
        class = "multischolar_absolute_path_in_artifact_state"
    )

    cycle <- oversized_bundle
    find_nested <- function(value) {
        if (identical(value$node_type, "nested_rectangular")) return(value)
        if (!identical(value$node_type, "list")) return(NULL)
        matches <- Filter(Negate(is.null), lapply(value$values, find_nested))
        if (length(matches) == 0L) NULL else matches[[1L]]
    }
    node <- find_nested(cycle$metadata$slot_values$args)
    expect_false(is.null(node))
    cycle_node <- list(
        node_type = "nested_rectangular",
        codec = list(
            id = .artifactNestedNodeCodec,
            version = .artifactNestedNodeCodecVersion
        ),
        payload_key = node$payload_key
    )
    encoded_cycle <- encodeArtifactTable(
        data.frame(
            serialized_node = as.character(jsonlite::serializeJSON(cycle_node)),
            stringsAsFactors = FALSE
        ),
        owner = "cycle"
    )
    cycle$payloads[[node$payload_key]] <- encoded_cycle$payload
    cycle$metadata$payloads[[node$payload_key]] <- encoded_cycle$metadata
    cycle$metadata$semantic_digest <- artifactSemanticDigest(
        artifactDiaBundleSemanticInput(cycle$metadata)
    )
    expect_error(
        hydrateDiaS4Artifact(cycle),
        class = "multischolar_invalid_artifact_payload"
    )
})

test_that("unknown codec versions and malformed payloads reject rather than guess", {
    bundle <- dehydrateDiaS4Artifact(makeDiaCodecObjects("iq")$peptide)
    future <- bundle
    future$metadata$schema_version <- 2L
    expect_error(
        hydrateDiaS4Artifact(future),
        class = "multischolar_unsupported_artifact_codec_version"
    )
    changed <- bundle
    changed$metadata$slot_values$args <- list(node_type = "null")
    expect_error(
        hydrateDiaS4Artifact(changed),
        class = "multischolar_artifact_semantic_digest_mismatch"
    )
    rectangular <- encodeArtifactTable(data.frame(id = 1:3))
    bad_shape <- rectangular$metadata
    bad_shape$dimensions$rows <- 4L
    expect_error(
        decodeArtifactRectangular(rectangular$payload, bad_shape),
        class = "multischolar_artifact_shape_mismatch"
    )
    future_table <- rectangular$metadata
    future_table$codec$version <- 2L
    expect_error(
        decodeArtifactRectangular(rectangular$payload, future_table),
        class = "multischolar_unsupported_artifact_codec_version"
    )
    wrong_schema <- rectangular$metadata
    wrong_schema$physical_schema[[1L]]$type <- "string"
    expect_error(
        decodeArtifactRectangular(rectangular$payload, wrong_schema),
        class = "multischolar_invalid_artifact_metadata"
    )
    unknown_field <- rectangular$metadata
    unknown_field$future_hint <- TRUE
    expect_error(
        decodeArtifactRectangular(rectangular$payload, unknown_field),
        class = "multischolar_invalid_artifact_metadata"
    )
})

test_that("codec helpers are collated once and introduce no production backend", {
    description <- read.dcf(test_path("..", "..", "DESCRIPTION"))
    collate <- strsplit(description[1L, "Collate"], "[[:space:]]+")[[1L]]
    collate <- gsub("^'|'$", "", collate[nzchar(collate)])
    helpers <- c(
        "utils_artifact_refs.R",
        "utils_artifact_table_codecs.R",
        "utils_artifact_codec_hydration.R",
        "utils_artifact_s4_codecs.R"
    )
    expect_true(all(vapply(helpers, function(helper) {
        sum(collate == helper) == 1L
    }, logical(1))))
    expect_lt(match("utils_artifact_paths.R", collate), match(helpers[[1L]], collate))
    expect_true(all(diff(match(helpers, collate)) > 0L))

    sources <- vapply(
        file.path(test_path("..", "..", "R"), helpers),
        function(path) paste(readLines(path, warn = FALSE), collapse = "\n"),
        character(1)
    )
    expect_false(any(grepl("DuckDB|dbConnect|reactiveVal|observeEvent", sources)))
    expect_false(any(grepl("saveRDS|readRDS", sources)))
})
