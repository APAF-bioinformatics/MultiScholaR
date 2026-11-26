# KIT-102 W2: the feature-table import core (R/utils_workflow_omics_import.R; the trx count policy
# and preview in R/func_trx_counts.R). Each moved function is checked against today's function on
# the module CI fixtures (MCI-012, MCI-018, MCI-026): equal payloads, refusal messages and artifact
# bytes. The parity blocks (text-numeric sample columns for lipid and trx) failed first
# (docs/qa/workflow-execution/kit/kit-102/w2/test-first.txt). The comparison blocks call today's
# `R/mod_*` import functions; when KIT-102 W5 deletes them, those blocks are quarantined or kept
# against frozen values (notes for phase C).

omicsImportConfig <- function(omic_type) {
    omicGuiConfig(omic_type)
}

omicsImportRequest <- function(assay, paths, assay_file, assay_name = "LCMS_Pos",
                               assay2_file = NULL, assay2_name = "", vendor_format,
                               detected_format = vendor_format,
                               id_col, annotation_col, sample_cols, sanitize_names = FALSE,
                               is_pattern = "^IS_") {
    list(assay1_name = assay_name, assay1_data = assay, assay1_file = assay_file,
        assay2_file = assay2_file, assay2_name = assay2_name, vendor_format = vendor_format,
        detected_format = detected_format, id_col = id_col, annotation_col = annotation_col,
        sample_columns = sample_cols, sanitize_names = sanitize_names, is_pattern = is_pattern,
        experiment_paths = paths)
}

omicsImportRun <- function(request, omic_type) {
    tryCatch(importOmicsTables(request, omicsImportConfig(omic_type)), error = identity)
}

# Paths under a run's project root, made relative, so two runs in two roots compare.
omicsImportRelative <- function(x, root) {
    rapply(x, function(value) gsub(root, "<root>", value, fixed = TRUE), classes = "character",
        how = "replace")
}

omicsImportLogWithoutTime <- function(log, root) {
    log$timestamp <- NULL
    omicsImportRelative(log, root)
}

# The import's artifacts (written paths and source copies), as relative path -> md5.
omicsImportArtifactBytes <- function(artifacts, root) {
    files <- c(unlist(artifacts$paths, use.names = FALSE), artifacts$source_copies)
    stats::setNames(unname(tools::md5sum(files)), gsub(root, "<root>", files, fixed = TRUE))
}

# KIT-116 W2 (Q-30): the comparison blocks compared the core with the pre-KIT-102 lipid, metab and trx
# import functions, which KIT-102 W5 deleted. Their old side is now frozen: literals, and sha256 digests
# of the canonical forms below, computed from the old functions at d469cea3 (kit-102/w2), where the
# comparison passed, by docs/qa/workflow-execution/kit/kit-116/freeze/q30-import.R (which also checked
# each old value against the core of that commit). A canonical form is the block's own comparison
# (expect_identical; paths under the run's root made relative). The core side stays live.
omicsImportDigest <- function(x) digest::digest(x, algo = "sha256")

# A successful import as the block compared it.
omicsImportRecord <- function(tables, sample_columns, data_format, column_mapping, log, validation,
                              artifacts, root) {
    list(tables = tables, sample_columns = sample_columns, data_format = data_format,
        column_mapping = column_mapping, log = omicsImportLogWithoutTime(log, root),
        validation = validation[c("valid", "errors", "warnings")],
        artifacts = omicsImportArtifactBytes(artifacts, root))
}

# A preview's fields as the block compared them, paths under `root` made relative.
omicsImportPreviewRecord <- function(headers, detected_format, format_confidence, import_result, id,
                                     annotation, is_pattern, assay_data, root) {
    omicsImportRelative(list(headers = headers, detected_format = detected_format,
        format_confidence = format_confidence, import_result = import_result, assay_data = assay_data,
        id = id, annotation = annotation, is_pattern = is_pattern), root)
}

# One case on the core against its frozen value: a refusal gives the frozen message and writes nothing;
# a success gives the frozen tables, mapping, log, validation and artifact bytes.
expectOmicsImportFrozen <- function(frozen, new, new_paths, info) {
    expect_false(is.null(frozen), info = info)
    if (!is.null(frozen$message)) {
        expect_true(inherits(new, "error"), info = paste(info, frozen$message))
        if (!inherits(new, "error")) return(invisible(NULL))
        expect_identical(conditionMessage(new), frozen$message, info = info)
        expect_length(list.files(new_paths$source_dir, recursive = TRUE), 0L)
        return(invisible(NULL))
    }
    expect_false(inherits(new, "error"), info = paste(info, if (inherits(new, "error")) {
        conditionMessage(new)
    }))
    if (inherits(new, "error")) return(invisible(NULL))
    expect_identical(omicsImportDigest(omicsImportRecord(new$assays, new$sample_columns,
        new$data_format, new$column_mapping, new$processing_log, new$validation, new$artifacts,
        new_paths$base_dir)), frozen$digest, info = info)
}

# The core's lipid import against the frozen pre-KIT-102 run (runLipidImportProcessing, through the
# MCI-018 helper).
compareLipidImport <- function(assay, info, assay_file = NULL, vendor_format = "lipidsearch",
                               detected_format = vendor_format, lipid_id_col = "LipidName",
                               annotation_col = "LipidClass",
                               sample_cols = module_ci_lipid_sample_cols(), assay_name = "LCMS_Pos",
                               assay2_file = NULL, assay2_name = "", sanitize_names = FALSE) {
    new_paths <- module_ci_lipid_import_paths()
    if (is.null(assay_file)) {
        assay_file <- module_ci_lipid_write_table(assay, file.path(new_paths$raw_dir,
            paste0(gsub("[/\\\\]", "_", assay_name), ".tsv")))
    }
    new <- omicsImportRun(omicsImportRequest(assay, new_paths, assay_file, assay_name, assay2_file,
        assay2_name, vendor_format, detected_format, lipid_id_col, annotation_col, sample_cols,
        sanitize_names), "lipidomics")
    expectOmicsImportFrozen(omicsImportFrozen_i5_lipid[[info]], new, new_paths, info)
    invisible(new)
}

# The core's metab import against the frozen pre-KIT-102 run (runMetabImportProcessing, through the
# MCI-012 helper).
compareMetabImport <- function(assay, info, assay_file = NULL, vendor_format = "custom",
                               detected_format = vendor_format, metabolite_col = "Feature.Name",
                               annotation_col = "Annotation",
                               sample_cols = grep("^(WT|KO|Sample|QC)", names(assay), value = TRUE),
                               assay_name = "LCMS_Pos", assay2_file = NULL, assay2_name = "",
                               sanitize_names = FALSE) {
    new_paths <- module_ci_metab_import_paths()
    if (is.null(assay_file)) {
        assay_file <- module_ci_metab_write_table(assay, file.path(new_paths$raw_dir,
            paste0(gsub("[/\\\\]", "_", assay_name), ".tsv")))
    }
    new <- omicsImportRun(omicsImportRequest(assay, new_paths, assay_file, assay_name, assay2_file,
        assay2_name, vendor_format, detected_format, metabolite_col, annotation_col, sample_cols,
        sanitize_names), "metabolomics")
    expectOmicsImportFrozen(omicsImportFrozen_i5_metab[[info]], new, new_paths, info)
    invisible(new)
}

# The core's trx import against the frozen pre-KIT-102 run (runTrxImportProcessing, through the MCI-026
# helper); the key is the fixture directory's name.
compareTrxImport <- function(counts_file, info) {
    new_paths <- module_ci_lipid_import_paths()
    preview <- previewTrxCountsImport(counts_file)
    new <- omicsImportRun(omicsImportRequest(preview$assay_data, new_paths, counts_file, "counts",
        NULL, "", "trx_counts", preview$detected_format, preview$updates$id$selected,
        preview$updates$annotation$selected, NULL, NULL, ""), "transcriptomics")
    expectOmicsImportFrozen(omicsImportFrozen_i5_trx[[info]], new, new_paths, info)
    invisible(new)
}

# ------------------------------------------------------------------------------------------------
# I1, I2: headers and preview.
# ------------------------------------------------------------------------------------------------

# The headers read by the pre-KIT-102 lipid (loadLipidImportAssayPreview) and metab
# (prepareMetabImportAssaySelectionState) default readers, frozen from d469cea3 (KIT-116).
omicsImportFrozen_i1_lipid <- list(
    tsv = c("Peak ID", "Name", "RT (min)", "Precursor m/z", "Adduct", "Formula", "Ontology", "Ion mode",
        "Total Score", "WT_1", "WT_2", "KO_1", "KO_2"),
    malformed = c("not", "a", "valid"),
    missing = character(0),
    csv = c("lipid_id", "annotation", "WT_1", "WT_2", "KO_1", "KO_2"))
omicsImportFrozen_i1_metab <- list(
    tsv = c("Peak ID", "Name", "RT (min)", "Precursor m/z", "Adduct", "Formula", "Ontology", "Ion mode",
        "Total Score", "WT_1", "WT_2", "KO_1", "KO_2"),
    malformed = c("not", "a", "valid"),
    missing = character(0),
    csv = c("lipid_id", "annotation", "WT_1", "WT_2", "KO_1", "KO_2"))

test_that("I1 the header reader is today's default reader (lipid and metab) (frozen from d469cea3, KIT-116)", {
    paths <- module_ci_lipid_import_paths()
    files <- c(
        tsv = module_ci_lipid_write_table(module_ci_lipid_msdial_assay(),
            file.path(paths$raw_dir, "a.tsv")),
        malformed = module_ci_lipid_write_malformed(file.path(paths$raw_dir, "m.tsv")),
        missing = file.path(paths$raw_dir, "absent.tsv"))
    csv <- file.path(paths$raw_dir, "q.csv")
    utils::write.csv(module_ci_lipid_custom_assay(), csv, row.names = FALSE)
    files <- c(files, csv = csv)
    for (name in names(files)) {
        headers <- suppressWarnings(readOmicsImportHeaders(files[[name]]))
        expect_identical(headers, omicsImportFrozen_i1_lipid[[name]], info = name)
        expect_identical(headers, omicsImportFrozen_i1_metab[[name]], info = name)
    }
    expect_identical(suppressWarnings(readOmicsImportHeaders(files[["missing"]])), character(0))
})

# The pre-KIT-102 previews (loadLipidImportAssayPreview, prepareMetabImportAssaySelectionState and
# trxImportPreview), as sha256 of omicsImportPreviewRecord, frozen from d469cea3 (KIT-116). The
# block compared import_result and the lipid assay_data with expect_equal; at d469cea3 they were
# identical, so the digest holds them to identity.
omicsImportFrozen_i2_lipid <- list(
    lipidsearch = "367c941b05e85813ade6b7da05bc3cc2bc37ab36af3dead72f35d328cccc1954",
    msdial = "d0c360c5b4195665e3209a6b1fa7b949d2d4007506276aeb018b4860ddba172c",
    no_annotation = "77962b1fa5f3d0ff30f46bc5e7128c255f8e959d453d0ef2884fd160cfd998d7",
    custom = "745d83b9d57b138838d5b75cee976b7d2aaa22bf399c39ae75dedc32faea1336",
    malformed = "01abb39641228ddb38702ca20e1d6a7ddc1f96dc7f31a700f65456cbd7a72857")
omicsImportFrozen_i2_metab <- list(
    msdial = "47a6b79b5821aa1dc72fc4e0040c40ef43b8227ce00bd01476bde4a630c42386",
    custom = "99fe68ea78469e229008a1538dd864ab34d9bfaea738a6208383ba98e112c54d")
omicsImportFrozen_i2_trx <- list(
    trx_screen = "4ddf739c1b81961996334b712d83a73281175be90d424227b5676c88683ebb53",
    counts_negative = "cc331cc3dac4fc7778710164b15c939c16b77262db33356a416a4a415b2c839b")

test_that("I2 the feature-table preview is today's lipid and metab preview (frozen from d469cea3, KIT-116)", {
    paths <- module_ci_lipid_import_paths()
    lipid_files <- c(
        lipidsearch = module_ci_lipid_write_table(module_ci_lipid_lipidsearch_assay(),
            file.path(paths$raw_dir, "ls.txt")),
        msdial = module_ci_lipid_write_table(module_ci_lipid_msdial_assay(),
            file.path(paths$raw_dir, "msdial_positive.tsv")),
        no_annotation = module_ci_lipid_write_table(
            module_ci_lipid_msdial_assay(include_annotation = FALSE),
            file.path(paths$raw_dir, "msdial_no_annotation.tsv")),
        custom = module_ci_lipid_write_table(module_ci_lipid_custom_assay(),
            file.path(paths$raw_dir, "custom.tsv")),
        malformed = module_ci_lipid_write_malformed(file.path(paths$raw_dir, "malformed.tsv")))
    lipid <- omicGuiConfig("lipidomics")$import
    for (name in names(lipid_files)) {
        new <- suppressWarnings(previewOmicsImport(lipid_files[[name]], lipid))
        expect_identical(omicsImportDigest(omicsImportPreviewRecord(new$headers, new$detected_format,
            new$format_confidence, new$import_result, new$updates$id, new$updates$annotation,
            new$updates$is_pattern, new$assay_data, paths$base_dir)), omicsImportFrozen_i2_lipid[[name]],
            info = name)
    }
    metab <- omicGuiConfig("metabolomics")$import
    metab_paths <- module_ci_metab_import_paths()
    metab_files <- c(
        msdial = module_ci_metab_write_table(module_ci_metab_msdial_assay(),
            file.path(metab_paths$raw_dir, "msdial_positive.tsv")),
        custom = module_ci_metab_write_table(module_ci_metab_custom_assay(),
            file.path(metab_paths$raw_dir, "custom.tsv")))
    for (name in names(metab_files)) {
        new <- suppressWarnings(previewOmicsImport(metab_files[[name]], metab))
        # The metab comparison did not include the assay data.
        expect_identical(omicsImportDigest(omicsImportPreviewRecord(new$headers, new$detected_format,
            new$format_confidence, new$import_result, new$updates$id, new$updates$annotation,
            new$updates$is_pattern, NULL, metab_paths$base_dir)), omicsImportFrozen_i2_metab[[name]],
            info = name)
    }
    absent <- file.path(paths$raw_dir, "absent.tsv")
    expect_error(suppressWarnings(previewOmicsImport(absent, lipid)),
        "Could not read headers from file", fixed = TRUE, class = "omics_import_invalid")
})

test_that("I2 the trx preview is today's trxImportPreview (frozen from d469cea3, KIT-116)", {
    for (dir in c(module_ci_trx_corpus(), module_ci_trx_fixture("counts_negative"))) {
        new <- previewTrxCountsImport(file.path(dir, "counts.tsv"))
        expect_identical(omicsImportDigest(omicsImportPreviewRecord(new$headers, new$detected_format,
            new$format_confidence, new$import_result, new$updates$id, new$updates$annotation, NULL,
            new$assay_data, dir)), omicsImportFrozen_i2_trx[[basename(dir)]], info = dir)
        expect_true("is_pattern" %in% names(new$updates))
        expect_null(new$updates$is_pattern)
    }
})

# ------------------------------------------------------------------------------------------------
# I3, I4: sample and column resolution.
# ------------------------------------------------------------------------------------------------

# The pre-KIT-102 resolutions, frozen from d469cea3 (KIT-116): resolveLipidImportSampleColumns (case,
# excludeNormalized), resolveMetabImportSampleColumns (case), resolveLipidImportEffectiveColumn and
# metab's column rule (resolveMetabImportColumnName for custom, else the selection) (vendor, name).
omicsImportFrozen_i3_lipid <- list(
    `1 FALSE` = c("Sample A", "Sample B", "Sample_normalized"),
    `1 TRUE` = c("Sample A", "Sample B", "Sample_normalized"),
    `2 FALSE` = c("Sample A", "Sample B", "QC-01", "QC-02", "Sample_normalized"),
    `2 TRUE` = c("Sample A", "Sample B", "QC-01", "QC-02"),
    `3 FALSE` = c("Sample A", "QC-01", "Sample_normalized"),
    `3 TRUE` = c("Sample A", "QC-01"),
    `4 FALSE` = c("Sample A", "QC-01", "Sample_normalized"),
    `4 TRUE` = c("Sample A", "QC-01"),
    `5 FALSE` = c("Sample A", "Sample B", "QC-01", "QC-02", "Sample_normalized"),
    `5 TRUE` = c("Sample A", "Sample B", "QC-01", "QC-02"))
omicsImportFrozen_i3_metab <- list(
    `1` = c("Sample A", "Sample B", "Sample_normalized"),
    `2` = c("Sample A", "Sample B", "QC-01", "QC-02", "Sample_normalized"),
    `3` = c("Sample A", "QC-01", "Sample_normalized"),
    `4` = c("Sample A", "QC-01", "Sample_normalized"),
    `5` = c("Sample A", "Sample B", "QC-01", "QC-02", "Sample_normalized"))
omicsImportFrozen_i4_lipid <- list(
    `custom LIPID_ID` = "lipid_id",
    `custom absent` = "absent",
    `custom ` = "",
    `custom NULL` = NULL,
    `msdial LIPID_ID` = "selected",
    `msdial absent` = "selected",
    `msdial ` = "selected",
    `msdial NULL` = "selected")
omicsImportFrozen_i4_metab <- list(
    `custom LIPID_ID` = "lipid_id",
    `custom absent` = "absent",
    `custom ` = "",
    `custom NULL` = NULL,
    `msdial LIPID_ID` = "selected",
    `msdial absent` = "selected",
    `msdial ` = "selected",
    `msdial NULL` = "selected")

test_that("I3, I4 the sample and column resolution is today's (lipid and metab) (frozen from d469cea3, KIT-116)", {
    custom <- module_ci_lipid_custom_assay(sample_prefix = "sample")
    custom$Sample_normalized <- 1
    import_result <- list(sample_columns = c("Sample A", "QC-01", "Sample_normalized"))
    cases <- list(
        list(vendor = "custom", pattern = "^sample", result = NULL),
        list(vendor = "custom", pattern = "^nomatch", result = NULL),
        list(vendor = "custom", pattern = "", result = import_result),
        list(vendor = "msdial", pattern = "^sample", result = import_result),
        list(vendor = "msdial", pattern = NULL, result = NULL))
    for (i in seq_along(cases)) {
        case <- cases[[i]]
        for (exclude in c(FALSE, TRUE)) {
            key <- paste(i, exclude)
            expect_true(key %in% names(omicsImportFrozen_i3_lipid), info = key)
            expect_identical(
                resolveOmicsImportSampleColumns(custom, case$vendor, case$pattern, case$result,
                    if (exclude) OMICS_IMPORT_NORMALISED_COLUMNS else NULL),
                omicsImportFrozen_i3_lipid[[key]], info = key)
        }
        expect_identical(
            resolveOmicsImportSampleColumns(custom, case$vendor, case$pattern, case$result),
            omicsImportFrozen_i3_metab[[as.character(i)]], info = i)
    }
    headers <- names(custom)
    for (vendor in c("custom", "msdial")) {
        for (custom_name in list("LIPID_ID", "absent", "", NULL)) {
            new <- resolveOmicsImportColumn(headers, vendor, "selected", custom_name)
            key <- paste(vendor, custom_name %||% "NULL")
            expect_true(key %in% names(omicsImportFrozen_i4_lipid), info = key)
            expect_identical(new, omicsImportFrozen_i4_lipid[[key]], info = key)
            expect_identical(new, omicsImportFrozen_i4_metab[[key]], info = key)
        }
    }
})

# ------------------------------------------------------------------------------------------------
# I5 to I13 through the orchestration, against today's runs.
# ------------------------------------------------------------------------------------------------

# The pre-KIT-102 lipid runs (runLipidImportProcessing through the MCI-018 helper), frozen from d469cea3
# (KIT-116): a refusal's message, or the sha256 of omicsImportRecord.
# KIT-118 M3 (2026-11-15): the digests of the cases whose fixtures hold an exact 0 in a
# sample column were re-frozen, because 0 is now read as missing with a warning: lipid
# `ls_zero`, `msdial`, `vendor_mismatch`, `pair`, `pair_unusual`; metabolomics `msdial`.
# Every other digest is unchanged.
omicsImportFrozen_i5_lipid <- list(
    ls_pos = list(digest = "b72fc8114f40fed8373d86ccce42746379304835fad0d0c3bb044f7cd4ca3a2e"),
    ls_zero = list(digest = "dffc369a5dba3f11dc541d20e3dd6c95fefd46a4079bebb077ad504a1072592f"),
    ls_no_class = list(digest = "fb6af95a0c1bfdc29374d1b4ff48ba31f778d617e5a70e729afd2f93d8b355e5"),
    ls_duplicate_ids = list(digest = "cb911979e5b0726a5ea1356e48f64f5050b4cc361d59a2865fcf699bf34af512"),
    missing_id = list(message = "Invalid lipidomics import: Assay 'LCMS_Pos': Lipid ID column not found: LipidName"),
    missing_id_column = list(
        message = "Invalid lipidomics import: Assay 'LCMS_Pos': Lipid ID column not found: MissingLipid"),
    no_samples = list(message = paste("Invalid lipidomics import: No sample columns specified;",
        "Assay 'LCMS_Pos': No sample columns specified")),
    non_numeric = list(message = "Invalid lipidomics import: Assay 'LCMS_Pos' has non-numeric sample columns: WT_1"),
    mismatched = list(
        message = "Invalid lipidomics import: Assay 'LCMS_Pos': Sample columns not found: Missing_1, Missing_2"),
    no_numeric = list(message = "Invalid lipidomics import: Assay 'LCMS_Pos' has non-numeric sample columns: WT_1"),
    msdial = list(digest = "e688ffde6df5b16e0fbef99911384e84e7d30d66d7cc8f9db23eac39daddcf3c"),
    vendor_mismatch = list(digest = "aa32ac9e94e4b9f63bb0773d38615a8841e3afe950aaa8f58879c21b44968070"),
    custom_override = list(digest = "db23fe88f5c761bf063125054b47fde2a643d6ff224399778dc1dd4087a96f19"),
    custom_null_detected = list(digest = "24a7c2966c6fa7492906ddc614d56a35ac8f33744e220c31b61522db0f8a2f2c"),
    path_name = list(message = "Assay names must not contain path separators"),
    empty_name = list(message = "Primary lipidomics assay name is required"),
    pair = list(digest = "26eb547810b2707b61f3252059f56ca13583284ff1012d320980e5aaac787d96"),
    pair_unusual = list(digest = "86efcb2344d992fa01ffd582d33c699ce60e2a771371175f1bd7cb6a0233e19c"),
    pair_duplicate = list(message = "Duplicate lipidomics assay names are not allowed"),
    pair_no_name = list(message = "Second lipidomics assay requires both a file and an assay name"),
    pair_no_file = list(message = "Second lipidomics assay requires both a file and an assay name"),
    pair_path = list(message = "Assay names must not contain path separators"),
    pair_empty_file = list(message = paste("Invalid lipidomics import: Assay 'LCMS_Neg': Lipid ID column not found:",
        "LipidName; Sample columns not found: WT_1, WT_2, KO_1, KO_2")),
    pair_msdial_reader = list(
        message = "Invalid lipidomics import: Assay 'LCMS_Neg': Lipid ID column not found: LipidName"))
# sanitizeLipidImportSampleNames: the sha256 of its assay list, and its sample columns.
omicsImportFrozen_i5_sanitize <- list(
    `TRUE` = list(assays = "697cb56fec881bdee42eda5705967b40a88e1a5cf120b877adeb461a963fda7b",
        sample_columns = c("sample_a", "sample_b", "qc_01", "qc_02")),
    `FALSE` = list(assays = "b246cd6e9108409d83bb767c905d915d0d527270217aec4405fe9ebc831b2fd5",
        sample_columns = c("Sample A", "Sample B", "QC-01", "QC-02")))

test_that("I5 to I13 the lipid import is today's on the MCI-018 matrix (frozen from d469cea3, KIT-116)", {
    skip_if_not(workflowOwnerCapable())
    paths <- module_ci_lipid_import_paths()
    for (case in list(
        list("ls_pos", module_ci_lipid_lipidsearch_assay()),
        list("ls_zero", module_ci_lipid_lipidsearch_assay(zero_values = TRUE),
            assay_name = "LCMS_Neg"),
        list("ls_no_class", module_ci_lipid_lipidsearch_assay(include_lipid_class = FALSE),
            annotation_col = ""),
        list("ls_duplicate_ids", module_ci_lipid_lipidsearch_assay(duplicate_lipid = TRUE)),
        list("missing_id", module_ci_lipid_lipidsearch_assay(include_lipid_name = FALSE)),
        list("missing_id_column", module_ci_lipid_lipidsearch_assay(),
            lipid_id_col = "MissingLipid"),
        list("no_samples", module_ci_lipid_lipidsearch_assay(), sample_cols = character()),
        list("non_numeric", module_ci_lipid_lipidsearch_assay(non_numeric_samples = TRUE)),
        list("mismatched", module_ci_lipid_lipidsearch_assay(),
            sample_cols = c("Missing_1", "Missing_2")),
        list("no_numeric", data.frame(LipidName = c("L1", "L2"), WT_1 = c("low", "high")),
            sample_cols = "WT_1"),
        list("msdial", module_ci_lipid_msdial_assay(duplicate_feature = TRUE),
            vendor_format = "msdial", lipid_id_col = "Peak ID", annotation_col = "Name"),
        list("vendor_mismatch", module_ci_lipid_msdial_assay(), vendor_format = "lipidsearch",
            detected_format = "msdial", lipid_id_col = "Peak ID", annotation_col = "Name"),
        list("custom_override", module_ci_lipid_lipidsearch_assay(), vendor_format = "custom",
            detected_format = "lipidsearch"),
        list("custom_null_detected", module_ci_lipid_custom_assay(), vendor_format = "custom",
            detected_format = NULL, lipid_id_col = "lipid_id", annotation_col = "annotation"),
        list("path_name", module_ci_lipid_lipidsearch_assay(), assay_name = "LCMS/Pos"),
        list("empty_name", module_ci_lipid_lipidsearch_assay(), assay_name = ""))) {
        args <- c(list(assay = case[[2]], info = case[[1]]), case[-(1:2)])
        do.call(compareLipidImport, args)
    }
    # Sanitising: today's lipid run notifies through shiny::showNotification, which the MCI helper's
    # stub session cannot show, so the step is compared directly (metab's runs end to end below).
    samples <- c("Sample A", "Sample B", "QC-01", "QC-02")
    custom <- module_ci_lipid_custom_assay(sample_prefix = "sample")
    assays <- list(a = custom, b = custom[c("lipid_id", "Sample A", "QC-02")])
    for (sanitize in c(TRUE, FALSE)) {
        frozen <- omicsImportFrozen_i5_sanitize[[as.character(sanitize)]]
        new <- sanitizeOmicsImportSampleNames(assays, samples, sanitize)
        expect_identical(omicsImportDigest(new$assays), frozen$assays, info = sanitize)
        expect_identical(new$sample_columns, frozen$sample_columns, info = sanitize)
        expect_identical(new$sanitized, sanitize)
    }
    # The second assay: read by the detected format's reader with the id and annotation columns.
    secondary <- module_ci_lipid_write_table(module_ci_lipid_lipidsearch_assay(zero_values = TRUE),
        file.path(paths$raw_dir, "lcms_neg.txt"))
    empty <- file.path(paths$raw_dir, "empty_assay2.txt")
    writeLines(character(), empty)
    primary <- module_ci_lipid_lipidsearch_assay()
    compareLipidImport(primary, "pair", assay2_file = secondary, assay2_name = "LCMS_Neg")
    compareLipidImport(primary, "pair_unusual", assay_name = "GCMS-2026_A", assay2_file = secondary,
        assay2_name = "LCMS_Pos_2026")
    compareLipidImport(primary, "pair_duplicate", assay2_file = secondary, assay2_name = "LCMS_Pos")
    compareLipidImport(primary, "pair_no_name", assay2_file = secondary, assay2_name = "")
    compareLipidImport(primary, "pair_no_file", assay2_file = NULL, assay2_name = "LCMS_Neg")
    compareLipidImport(primary, "pair_path", assay2_file = secondary, assay2_name = "a\\b")
    compareLipidImport(primary, "pair_empty_file", assay2_file = empty, assay2_name = "LCMS_Neg")
    msdial_neg <- module_ci_lipid_write_table(module_ci_lipid_msdial_assay(),
        file.path(paths$raw_dir, "msdial_neg.tsv"))
    compareLipidImport(primary, "pair_msdial_reader", vendor_format = "msdial",
        detected_format = "msdial", assay2_file = msdial_neg, assay2_name = "LCMS_Neg")
})

# The pre-KIT-102 metab runs (runMetabImportProcessing through the MCI-012 helper), frozen from
# d469cea3 (KIT-116): a refusal's message, or the sha256 of omicsImportRecord.
omicsImportFrozen_i5_metab <- list(
    annotation_present = list(digest = "fb16f3b5fa21456363a6b177d521427711c65b8acd97621ee16eec1be14b9e68"),
    annotation_absent = list(digest = "efe3283a6676cf7a2485b8f0cfec0fe2e8b8cf43fb66333bb050c0a630da3e55"),
    text_numeric = list(digest = "ef7034f51eff11297e19f50b766957158eb74cddf46152c148d49b0359cf44b6"),
    msdial = list(digest = "f244efc758d4e01d113eb839b37802f2ef452cbac1a8995947fb0a7db7c87850"),
    missing_id_column = list(message = "Invalid metabolomics import: Metabolite ID column not found: Missing.ID"),
    empty_id = list(message = paste("Invalid metabolomics import: Metabolite ID column is not specified.;",
        "Metabolite ID column is not specified")),
    no_samples = list(message = "Invalid metabolomics import: No sample columns specified.; No sample columns specified"),
    non_numeric = list(
        message = "Invalid metabolomics import: Assay 'LCMS_Pos' has non-numeric sample column(s): WT_1"),
    infinite = list(message = "Invalid metabolomics import: Assay 'LCMS_Pos' has non-numeric sample column(s): WT_1"),
    mismatched = list(message = "Invalid metabolomics import: Sample columns not found: Missing_1"),
    malformed = list(message = paste("Invalid metabolomics import: Metabolite ID column not found: Feature.Name;",
        "Sample columns not found: WT_1")),
    sanitised = list(digest = "a00f11fdace064662eb10e18ab0a31edd79c65f94e18fc4414324ea6df25878e"),
    path_name = list(message = paste("Invalid metabolomics import: Assay names must not contain path separators",
        "or control characters.")),
    control_name = list(message = paste("Invalid metabolomics import: Assay names must not contain path separators",
        "or control characters.")),
    empty_name = list(message = "Metabolomics import requires non-empty assay names."),
    pair = list(digest = "83c534c6ad93db05247d3620c66a6fe479eab059d12c10ce6e88522846536215"),
    pair_unusual = list(digest = "9a2e963b7abcd989dfa3dd947e8e64aee8667eb811f4c85144929e016d40fb9e"),
    pair_duplicate = list(message = "Duplicate metabolomics assay names are not allowed: LCMS_Pos"),
    pair_no_name = list(digest = "331658c189fa97ae57cf3f41fb2047f83aa1f266490f154cfb8a35de478ed5b6"),
    pair_path = list(message = paste("Invalid metabolomics import: Assay names must not contain path separators",
        "or control characters.")))

test_that("I5 to I13 the metab import is today's on the MCI-012 matrix (frozen from d469cea3, KIT-116)", {
    skip_if_not(workflowOwnerCapable())
    paths <- module_ci_metab_import_paths()
    valid <- module_ci_metab_custom_assay()
    for (case in list(
        list("annotation_present", module_ci_metab_custom_assay(numeric_as_character = FALSE)),
        list("annotation_absent", module_ci_metab_custom_assay(annotation = FALSE),
            annotation_col = ""),
        list("text_numeric",
            module_ci_metab_custom_assay(numeric_as_character = TRUE, sample_prefix = "KO"),
            sample_cols = c("Sample-A", "Sample-B", "QC-01", "QC-02")),
        list("msdial", module_ci_metab_msdial_assay(duplicate_feature = TRUE),
            vendor_format = "msdial", metabolite_col = "Peak ID", annotation_col = "Name",
            sample_cols = c("WT_1", "WT_2", "KO_1", "KO_2")),
        list("missing_id_column", valid, metabolite_col = "Missing.ID",
            sample_cols = c("WT_1", "WT_2")),
        list("empty_id", valid, metabolite_col = "", sample_cols = c("WT_1", "WT_2")),
        list("no_samples", valid, sample_cols = character()),
        list("non_numeric", within(valid, WT_1 <- c("ok", "bad", "values")),
            sample_cols = c("WT_1", "WT_2")),
        list("infinite", within(valid, WT_1 <- c(1, Inf, 3)), sample_cols = c("WT_1", "WT_2")),
        list("mismatched", valid, sample_cols = "Missing_1"),
        list("malformed", data.frame(Bad = c("{", "}")), sample_cols = "WT_1"),
        list("sanitised", module_ci_metab_custom_assay(sample_prefix = "KO"),
            sample_cols = c("Sample-A", "Sample-B", "QC-01", "QC-02"), sanitize_names = TRUE),
        list("path_name", valid, assay_name = "LCMS/Pos"),
        list("control_name", valid, assay_name = "LCMS\tPos"),
        list("empty_name", valid, assay_name = ""))) {
        args <- c(list(assay = case[[2]], info = case[[1]]), case[-(1:2)])
        do.call(compareMetabImport, args)
    }
    secondary <- module_ci_metab_write_table(module_ci_metab_custom_assay(),
        file.path(paths$raw_dir, "lcms_neg.tsv"))
    compareMetabImport(valid, "pair", assay2_file = secondary, assay2_name = "LCMS_Neg")
    compareMetabImport(valid, "pair_unusual", assay_name = "GCMS-2026_A", assay2_file = secondary,
        assay2_name = "LCMS_Pos_2026")
    compareMetabImport(valid, "pair_duplicate", assay2_file = secondary, assay2_name = "LCMS_Pos")
    compareMetabImport(valid, "pair_no_name", assay_name = "GCMS", assay2_file = secondary,
        assay2_name = "")
    compareMetabImport(valid, "pair_path", assay2_file = secondary, assay2_name = "LCMS\nNeg")
})

# The pre-KIT-102 trx runs (runTrxImportProcessing through the MCI-026 helper), frozen from d469cea3
# (KIT-116): a refusal's message, or the sha256 of omicsImportRecord.
omicsImportFrozen_i5_trx <- list(
    trx_screen = list(digest = "5f44b8b14a194af5e365b1722a325767545d1249e18a13625eed7f6890b6e04d"),
    counts_whole_double = list(digest = "c694fe5c0e2286ce7e2fdc0edbbcc436c61d6719495026d1340871f3192cf265"),
    counts_na = list(message = "NA or NaN count (1 cell(s); first: gene 1, sample RPMI_03_r2) (trx_counts_na)"),
    counts_negative = list(
        message = "negative count (1 cell(s); first: gene 1, sample RPMI_03_r2) (trx_counts_negative)"),
    counts_fractional = list(
        message = "non-whole count (1 cell(s); first: gene 1, sample RPMI_03_r2) (trx_counts_fractional)"))

test_that("I5 to I13 the trx import is today's on the MCI-026 corpus and fixtures (frozen from d469cea3, KIT-116)", {
    skip_if_not(workflowOwnerCapable())
    for (dir in c(module_ci_trx_corpus(), module_ci_trx_fixture("counts_whole_double"),
        module_ci_trx_fixture("counts_na"), module_ci_trx_fixture("counts_negative"),
        module_ci_trx_fixture("counts_fractional"))) {
        compareTrxImport(file.path(dir, "counts.tsv"), basename(dir))
    }
    # Refused by the count policy's class, as today.
    refused <- c(counts_na = "trx_counts_na", counts_negative = "trx_counts_negative",
        counts_fractional = "trx_counts_fractional")
    for (name in names(refused)) {
        file <- file.path(module_ci_trx_fixture(name), "counts.tsv")
        preview <- previewTrxCountsImport(file)
        request <- omicsImportRequest(preview$assay_data, module_ci_lipid_import_paths(), file,
            "counts", vendor_format = "trx_counts", id_col = "GeneID",
            annotation_col = "GeneName", sample_cols = NULL)
        expect_error(importOmicsTables(request, omicsImportConfig("transcriptomics")),
            class = refused[[name]])
    }
})

test_that("the import adds no copy of an assay that needs no conversion (scoping risk 4)", {
    assay <- module_ci_lipid_lipidsearch_assay()
    paths <- module_ci_lipid_import_paths()
    request <- omicsImportRequest(assay, paths, NULL, vendor_format = "lipidsearch",
        id_col = "LipidName", annotation_col = "LipidClass",
        sample_cols = module_ci_lipid_sample_cols())
    payload <- importOmicsTables(request, omicsImportConfig("lipidomics"))
    expect_identical(rlang::obj_address(payload$assays$LCMS_Pos), rlang::obj_address(assay))
})

# ------------------------------------------------------------------------------------------------
# The shared conversion helper and the parity (failed first).
# ------------------------------------------------------------------------------------------------

test_that("coerceOmicsNumericLikeColumns converts only numeric-like text and factor columns (frozen from d469cea3, KIT-116)", {
    table <- data.frame(id = c("a", "b", "c"), text = c(" 1.5", "", NA),
        factor = factor(c("2", "3", "4")), word = c("1", "x", "2"), number = c(1, 2, 3),
        flag = c(NA, NA, NA), unlisted = c("7", "8", "9"), infinite = c("1", "Inf", "2"),
        stringsAsFactors = FALSE)
    out <- coerceOmicsNumericLikeColumns(table, c("text", "factor", "word", "number", "flag",
        "infinite", "absent"))
    expect_identical(out$text, c(1.5, NA, NA))
    expect_identical(out$factor, c(2, 3, 4))
    expect_identical(out$word, table$word)
    expect_identical(out$number, table$number)
    expect_identical(out$flag, table$flag)
    expect_identical(out$unlisted, table$unlisted)
    expect_identical(out$infinite, table$infinite)
    expect_identical(out$id, table$id)
    # Nothing to convert: the same object.
    numeric_only <- table[c("id", "number")]
    expect_identical(rlang::obj_address(coerceOmicsNumericLikeColumns(numeric_only, "number")),
        rlang::obj_address(numeric_only))
    # The rule is metab's (metabImportValuesAreNumericLike, coerceMetabImportWorkflowPayloadSamples);
    # their results on this table, frozen from d469cea3 (KIT-116).
    numeric_like <- c(id = FALSE, text = TRUE, factor = TRUE, word = FALSE, number = TRUE, flag = FALSE,
        unlisted = TRUE, infinite = FALSE)
    for (column in names(table)) {
        expect_identical(omicsImportValuesAreNumericLike(table[[column]]), numeric_like[[column]],
            info = column)
    }
    small <- table[c("id", "text", "factor")]
    payload <- structure(list(id = c("a", "b", "c"), text = c(1.5, NA, NA), factor = c(2, 3, 4)),
        row.names = c(NA, -3L), class = "data.frame")
    expect_identical(coerceOmicsNumericLikeColumns(small, c("text", "factor")), payload)
})

test_that("parity: lipid text-numeric sample columns are converted, not refused (failed first)", {
    paths <- module_ci_lipid_import_paths()
    text <- module_ci_lipid_custom_assay(numeric_as_character = TRUE)
    number <- module_ci_lipid_custom_assay(numeric_as_character = FALSE)
    request <- function(assay, root) {
        omicsImportRequest(assay, root, NULL, vendor_format = "custom", detected_format = NULL,
            id_col = "lipid_id", annotation_col = "annotation",
            sample_cols = module_ci_lipid_sample_cols())
    }
    config <- omicsImportConfig("lipidomics")
    converted <- importOmicsTables(request(text, paths), config)
    samples <- converted$assays$LCMS_Pos[module_ci_lipid_sample_cols()]
    expect_true(all(vapply(samples, is.double, logical(1))))
    from_numbers <- importOmicsTables(request(number, module_ci_lipid_import_paths()), config)
    expect_identical(converted$assays, from_numbers$assays)
    # A column that is not numeric-like is still refused, with today's message.
    text$WT_1 <- c("ok", "bad", "values")
    expect_error(importOmicsTables(request(text, module_ci_lipid_import_paths()), config),
        "Invalid lipidomics import: Assay 'LCMS_Pos' has non-numeric sample columns: WT_1",
        fixed = TRUE, class = "omics_import_invalid")
})

test_that("parity: trx converts text-numeric counts before its count policy (failed first)", {
    config <- omicsImportConfig("transcriptomics")
    run <- function(values) {
        table <- data.frame(GeneID = c("G1", "G2"), GeneName = c("N1", "N2"), S1 = values,
            S2 = c(5, 7), stringsAsFactors = FALSE)
        importOmicsTables(omicsImportRequest(table, module_ci_lipid_import_paths(), NULL, "counts",
            vendor_format = "trx_counts", id_col = "GeneID", annotation_col = "GeneName",
            sample_cols = NULL, is_pattern = ""), config)
    }
    passed <- run(c("12", "3"))
    expect_identical(passed$assays$counts$S1, c(12, 3))
    expect_identical(passed$sample_columns, c("S1", "S2"))
    expect_error(run(c("-1", "3")), "negative count", class = "trx_counts_negative")
    expect_error(run(c("x", "3")), "counts must be numeric (trx_counts_invalid)", fixed = TRUE,
        class = "trx_counts_invalid")
})

test_that("checkTrxImportCounts keeps its name and refusals in the trx core (MCI-026.2)", {
    table <- utils::read.delim(file.path(module_ci_trx_corpus(), "counts.tsv"), check.names = FALSE,
        stringsAsFactors = FALSE)
    expect_identical(checkTrxImportCounts(list(counts = table), "GeneID", "GeneName"),
        list(counts = table))
    expect_error(checkTrxImportCounts(list(counts = table), "GeneID", ""),
        class = "trx_counts_invalid")
    expect_error(checkTrxImportCounts(list(counts = table), "GeneID", "Absent"),
        "lacks the column(s) Absent", fixed = TRUE, class = "trx_counts_invalid")
    table$note <- "x"
    expect_error(checkTrxImportCounts(list(counts = table), "GeneID", "GeneName"),
        class = "trx_counts_invalid")
})

test_that("F102-01: metab's missing-value warning counts empty text cells, after the conversion", {
    # Today's metab validation ran before the conversion, so an empty text cell was neither NA nor
    # 0 and the warning missed it (test-first.txt). The tables are unchanged.
    assay <- module_ci_metab_custom_assay(numeric_as_character = TRUE)
    for (column in c("WT_1", "WT_2")) assay[[column]] <- c("", "", "")
    assay$KO_1[1] <- ""
    payload <- importOmicsTables(omicsImportRequest(assay, module_ci_metab_import_paths(), NULL,
        vendor_format = "custom", id_col = "Feature.Name", annotation_col = "Annotation",
        sample_cols = c("WT_1", "WT_2", "KO_1", "KO_2")), omicsImportConfig("metabolomics"))
    expect_true(payload$validation$valid)
    expect_identical(payload$validation$warnings, "High proportion of missing values: 58.3%")
    expect_identical(payload$assays$LCMS_Pos$WT_1, rep(NA_real_, 3))
})

# ------------------------------------------------------------------------------------------------
# KIT-118 M3 (KIT-114 (e); user's call 2026-11-15): for lipid and metab a 0 in a sample column is "not detected", read as
# missing with its count reported; a negative refuses the import. Transcriptomics counts keep their zeros.
# ------------------------------------------------------------------------------------------------

test_that("(KIT-118 M3) lipid and metab read a zero as missing and report the count; a negative refuses (fails before)", {
    samples <- module_ci_lipid_sample_cols()
    request <- function(assay) omicsImportRequest(assay, module_ci_lipid_import_paths(), NULL, vendor_format = "custom",
        detected_format = NULL, id_col = "lipid_id", annotation_col = "annotation", sample_cols = samples)
    for (omic in c("lipidomics", "metabolomics")) {
        config <- omicsImportConfig(omic)
        clean <- module_ci_lipid_custom_assay(numeric_as_character = FALSE)
        untouched <- importOmicsTables(request(clean), config)
        expect_identical(untouched$assays[[1L]][samples], clean[samples], info = omic)
        expect_null(untouched$zeros_as_missing, info = omic)
        with_zero <- clean
        with_zero[[samples[[1L]]]][c(1L, 2L)] <- 0
        with_zero[[samples[[2L]]]][[3L]] <- 0
        payload <- importOmicsTables(request(with_zero), config)
        got <- payload$assays[[1L]]
        expect_identical(which(is.na(got[[samples[[1L]]]])), c(1L, 2L), info = omic)
        expect_identical(which(is.na(got[[samples[[2L]]]])), 3L, info = omic)
        expected <- clean[samples]
        expected[[samples[[1L]]]][c(1L, 2L)] <- NA_real_
        expected[[samples[[2L]]]][[3L]] <- NA_real_
        expect_identical(got[samples], expected, info = omic)
        expect_identical(payload$zeros_as_missing, stats::setNames(3L, names(payload$assays)[[1L]]), info = omic)
        expect_true(any(grepl("3 zero value(s) read as missing", payload$validation$warnings, fixed = TRUE)), info = omic)
        negative <- clean
        negative[[samples[[1L]]]][[1L]] <- -5
        expect_error(importOmicsTables(request(negative), config), "1 negative value(s)", fixed = TRUE,
            class = "omics_import_invalid")
    }
})

test_that("(KIT-118 M3) transcriptomics counts keep their zeros", {
    table <- data.frame(GeneID = c("G1", "G2", "G3"), GeneName = c("N1", "N2", "N3"), S1 = c(0, 4, 9), S2 = c(5, 0, 7),
        stringsAsFactors = FALSE)
    payload <- importOmicsTables(omicsImportRequest(table, module_ci_lipid_import_paths(), NULL, "counts",
        vendor_format = "trx_counts", id_col = "GeneID", annotation_col = "GeneName", sample_cols = NULL, is_pattern = ""),
        omicsImportConfig("transcriptomics"))
    expect_identical(payload$assays$counts$S1, c(0, 4, 9))
    expect_null(payload$zeros_as_missing)
})
