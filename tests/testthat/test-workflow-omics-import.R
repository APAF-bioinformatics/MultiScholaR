# KIT-102 W2: the feature-table import core (R/utils_workflow_omics_import.R; the trx count policy
# and preview in R/func_trx_counts.R). Each moved function is checked against today's function on
# the module CI fixtures (MCI-012, MCI-018, MCI-026): equal payloads, refusal messages and artifact
# bytes. The parity blocks (text-numeric sample columns for lipid and trx) failed first
# (docs/qa/workflow-execution/kit/kit-102/w2/test-first.txt). The comparison blocks call today's
# `R/mod_*` import functions; when KIT-102 W5 deletes them, those blocks are quarantined or kept
# against frozen values (notes for phase C).

omicsImportConfig <- function(omic_type) {
    config <- omicGuiConfig(omic_type)
    # The lead wires the trx `check` plugin into the configuration (W1 notes); set it here until
    # then.
    if (identical(omic_type, "transcriptomics")) config$import$check <- checkTrxImportCounts
    config
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

# One case on today's path and on the core: a refusal gives the same message and writes nothing; a
# success gives the same tables, mapping, log, validation and artifact bytes.
expectOmicsImportLikeToday <- function(old, new, old_root, new_paths, info) {
    if (is.null(old$tables)) {
        expect_true(inherits(new, "error"), info = paste(info, old$message))
        if (!inherits(new, "error")) return(invisible(NULL))
        expect_identical(conditionMessage(new), old$message, info = info)
        expect_length(list.files(new_paths$source_dir, recursive = TRUE), 0L)
        return(invisible(NULL))
    }
    expect_false(inherits(new, "error"), info = paste(info, if (inherits(new, "error")) {
        conditionMessage(new)
    }))
    if (inherits(new, "error")) return(invisible(NULL))
    expect_identical(new$assays, old$tables, info = info)
    expect_identical(new$sample_columns, old$sample_columns, info = info)
    expect_identical(new$data_format, old$data_format, info = info)
    expect_identical(new$column_mapping, old$column_mapping, info = info)
    expect_identical(omicsImportLogWithoutTime(new$processing_log, new_paths$base_dir),
        omicsImportLogWithoutTime(old$log, old_root), info = info)
    expect_identical(new$validation[c("valid", "errors", "warnings")],
        old$validation[c("valid", "errors", "warnings")], info = info)
    expect_identical(omicsImportArtifactBytes(new$artifacts, new_paths$base_dir),
        omicsImportArtifactBytes(old$artifacts, old_root), info = info)
}

# Today's lipid import (runLipidImportProcessing, through the MCI-018 helper) and the core.
compareLipidImport <- function(assay, info, assay_file = NULL, vendor_format = "lipidsearch",
                               detected_format = vendor_format, lipid_id_col = "LipidName",
                               annotation_col = "LipidClass",
                               sample_cols = module_ci_lipid_sample_cols(), assay_name = "LCMS_Pos",
                               assay2_file = NULL, assay2_name = "", sanitize_names = FALSE) {
    old_paths <- module_ci_lipid_import_paths()
    if (is.null(assay_file)) {
        assay_file <- module_ci_lipid_write_table(assay, file.path(old_paths$raw_dir,
            paste0(gsub("[/\\\\]", "_", assay_name), ".tsv")))
    }
    old <- module_ci_lipid_import_run(assay = assay, paths = old_paths, assay_name = assay_name,
        assay_file = assay_file, assay2_file = assay2_file, assay2_name = assay2_name,
        vendor_format = vendor_format, detected_format = detected_format,
        lipid_id_col = lipid_id_col, annotation_col = annotation_col, sample_cols = sample_cols,
        sanitize_names = sanitize_names)
    today <- if (is.null(old$result)) {
        list(message = sub("^Error processing import: ", "", old$error_messages[[1]]))
    } else {
        list(tables = old$result$assayList, sample_columns = old$result$sampleColumns,
            data_format = old$workflow$data_format, column_mapping = old$workflow$column_mapping,
            log = old$workflow$processing_log$setup_import,
            validation = old$result$validationResult, artifacts = old$result$artifactResult)
    }
    new_paths <- module_ci_lipid_import_paths()
    new <- omicsImportRun(omicsImportRequest(assay, new_paths, assay_file, assay_name, assay2_file,
        assay2_name, vendor_format, detected_format, lipid_id_col, annotation_col, sample_cols,
        sanitize_names), "lipidomics")
    expectOmicsImportLikeToday(today, new, old_paths$base_dir, new_paths, info)
    invisible(new)
}

# Today's metab import (runMetabImportProcessing, through the MCI-012 helper) and the core.
compareMetabImport <- function(assay, info, assay_file = NULL, vendor_format = "custom",
                               detected_format = vendor_format, metabolite_col = "Feature.Name",
                               annotation_col = "Annotation",
                               sample_cols = grep("^(WT|KO|Sample|QC)", names(assay), value = TRUE),
                               assay_name = "LCMS_Pos", assay2_file = NULL, assay2_name = "",
                               sanitize_names = FALSE) {
    old_paths <- module_ci_metab_import_paths()
    if (is.null(assay_file)) {
        assay_file <- module_ci_metab_write_table(assay, file.path(old_paths$raw_dir,
            paste0(gsub("[/\\\\]", "_", assay_name), ".tsv")))
    }
    old <- module_ci_metab_import_run(assay = assay, paths = old_paths, assay_name = assay_name,
        assay_file = assay_file, assay2_file = assay2_file, assay2_name = assay2_name,
        vendor_format = vendor_format, detected_format = detected_format,
        metabolite_col = metabolite_col, annotation_col = annotation_col,
        sample_cols = sample_cols, sanitize_names = sanitize_names)
    today <- if (identical(old$result$status, "error")) {
        list(message = conditionMessage(old$result$error))
    } else {
        payload <- old$result$workflowPayload
        list(tables = payload$assayList, sample_columns = payload$sampleCols,
            data_format = old$workflow$data_format, column_mapping = old$workflow$column_mapping,
            log = old$workflow$processing_log$setup_import,
            validation = old$result$validationResult, artifacts = old$result$artifactResult)
    }
    new_paths <- module_ci_metab_import_paths()
    new <- omicsImportRun(omicsImportRequest(assay, new_paths, assay_file, assay_name, assay2_file,
        assay2_name, vendor_format, detected_format, metabolite_col, annotation_col, sample_cols,
        sanitize_names), "metabolomics")
    expectOmicsImportLikeToday(today, new, old_paths$base_dir, new_paths, info)
    invisible(new)
}

# Today's trx import (runTrxImportProcessing, through the MCI-026 helper) and the core.
compareTrxImport <- function(counts_file, info) {
    fixture <- module_ci_trx_fixture_owner()
    old <- module_ci_trx_import(fixture, counts_file)
    today <- if (is.null(old$result)) {
        list(message = sub("^Error processing import: ", "", old$errors[[1]]))
    } else {
        list(tables = old$result$assayList, sample_columns = old$result$sampleColumns,
            data_format = fixture$workflow$data_format,
            column_mapping = fixture$workflow$column_mapping,
            log = fixture$workflow$processing_log$setup_import,
            validation = old$result$validationResult, artifacts = old$result$artifactResult)
    }
    new_paths <- module_ci_lipid_import_paths()
    preview <- previewTrxCountsImport(counts_file)
    new <- omicsImportRun(omicsImportRequest(preview$assay_data, new_paths, counts_file, "counts",
        NULL, "", "trx_counts", preview$detected_format, preview$updates$id$selected,
        preview$updates$annotation$selected, NULL, NULL, ""), "transcriptomics")
    expectOmicsImportLikeToday(today, new, fixture$paths$base_dir, new_paths, info)
    invisible(new)
}

# ------------------------------------------------------------------------------------------------
# I1, I2: headers and preview.
# ------------------------------------------------------------------------------------------------

test_that("I1 the header reader is today's default reader (lipid and metab)", {
    lipid_default <- eval(formals(loadLipidImportAssayPreview)$readHeaders)
    metab_default <- eval(formals(prepareMetabImportAssaySelectionState)$readHeadersFn)
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
        expect_identical(headers, suppressWarnings(lipid_default(files[[name]])), info = name)
        expect_identical(headers, suppressWarnings(metab_default(files[[name]])), info = name)
    }
    expect_identical(suppressWarnings(readOmicsImportHeaders(files[["missing"]])), character(0))
})

test_that("I2 the feature-table preview is today's lipid and metab preview", {
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
        old <- suppressWarnings(loadLipidImportAssayPreview(lipid_files[[name]],
            logInfo = function(...) NULL))
        new <- suppressWarnings(previewOmicsImport(lipid_files[[name]], lipid))
        expect_identical(new$headers, old$headers, info = name)
        expect_identical(new$detected_format, old$detectedFormat, info = name)
        expect_identical(new$format_confidence, old$formatConfidence, info = name)
        expect_equal(new$import_result, old$importResult, info = name)
        expect_equal(new$assay_data, old$assayData, info = name)
        expect_identical(new$updates$id, old$updates$lipidId, info = name)
        expect_identical(new$updates$annotation, old$updates$annotation, info = name)
        expect_identical(new$updates$is_pattern, old$updates$isPattern, info = name)
    }
    metab <- omicGuiConfig("metabolomics")$import
    metab_paths <- module_ci_metab_import_paths()
    metab_files <- c(
        msdial = module_ci_metab_write_table(module_ci_metab_msdial_assay(),
            file.path(metab_paths$raw_dir, "msdial_positive.tsv")),
        custom = module_ci_metab_write_table(module_ci_metab_custom_assay(),
            file.path(metab_paths$raw_dir, "custom.tsv")))
    for (name in names(metab_files)) {
        old <- suppressWarnings(prepareMetabImportAssaySelectionState(metab_files[[name]]))
        new <- suppressWarnings(previewOmicsImport(metab_files[[name]], metab))
        expect_identical(new$headers, old$headers, info = name)
        expect_identical(list(format = new$detected_format, confidence = new$format_confidence),
            old$formatInfo[c("format", "confidence")], info = name)
        expect_equal(new$import_result, old$importResult, info = name)
        expect_identical(new$updates$id, list(choices = old$metaboliteIdChoices,
            selected = old$selectedMetaboliteId), info = name)
        expect_identical(new$updates$annotation, list(choices = old$annotationChoices,
            selected = old$selectedAnnotation), info = name)
        expect_identical(new$updates$is_pattern, old$isPattern, info = name)
    }
    absent <- file.path(paths$raw_dir, "absent.tsv")
    expect_error(suppressWarnings(previewOmicsImport(absent, lipid)),
        "Could not read headers from file", fixed = TRUE, class = "omics_import_invalid")
})

test_that("I2 the trx preview is today's trxImportPreview", {
    for (dir in c(module_ci_trx_corpus(), module_ci_trx_fixture("counts_negative"))) {
        old <- trxImportPreview(file.path(dir, "counts.tsv"))
        new <- previewTrxCountsImport(file.path(dir, "counts.tsv"))
        expect_identical(new$headers, old$headers, info = dir)
        expect_identical(new$detected_format, old$detectedFormat, info = dir)
        expect_identical(new$format_confidence, old$formatConfidence, info = dir)
        expect_identical(new$import_result, old$importResult, info = dir)
        expect_identical(new$assay_data, old$assayData, info = dir)
        expect_identical(new$updates$id, old$updates$lipidId, info = dir)
        expect_identical(new$updates$annotation, old$updates$annotation, info = dir)
        expect_true("is_pattern" %in% names(new$updates))
        expect_null(new$updates$is_pattern)
    }
})

# ------------------------------------------------------------------------------------------------
# I3, I4: sample and column resolution.
# ------------------------------------------------------------------------------------------------

test_that("I3, I4 the sample and column resolution is today's (lipid and metab)", {
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
            expect_identical(
                resolveOmicsImportSampleColumns(custom, case$vendor, case$pattern, case$result,
                    if (exclude) OMICS_IMPORT_NORMALISED_COLUMNS else NULL),
                resolveLipidImportSampleColumns(custom, case$result, case$vendor,
                    case$pattern %||% "", excludeNormalized = exclude), info = paste(i, exclude))
        }
        expect_identical(
            resolveOmicsImportSampleColumns(custom, case$vendor, case$pattern, case$result),
            resolveMetabImportSampleColumns(custom, case$vendor, case$pattern, case$result),
            info = i)
    }
    headers <- names(custom)
    for (vendor in c("custom", "msdial")) {
        for (custom_name in list("LIPID_ID", "absent", "", NULL)) {
            new <- resolveOmicsImportColumn(headers, vendor, "selected", custom_name)
            expect_identical(new, resolveLipidImportEffectiveColumn(custom, vendor, "selected",
                custom_name), info = paste(vendor, custom_name %||% "NULL"))
            metab_today <- if (identical(vendor, "custom")) {
                resolveMetabImportColumnName(headers, custom_name)
            } else {
                "selected"
            }
            expect_identical(new, metab_today, info = paste(vendor, custom_name %||% "NULL"))
        }
    }
})

# ------------------------------------------------------------------------------------------------
# I5 to I13 through the orchestration, against today's runs.
# ------------------------------------------------------------------------------------------------

test_that("I5 to I13 the lipid import is today's on the MCI-018 matrix", {
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
        old <- sanitizeLipidImportSampleNames(assays, samples, sanitize,
            logInfo = function(...) NULL, notify = function(...) NULL)
        new <- sanitizeOmicsImportSampleNames(assays, samples, sanitize)
        expect_identical(new$assays, old$assayList, info = sanitize)
        expect_identical(new$sample_columns, old$sampleColumns, info = sanitize)
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

test_that("I5 to I13 the metab import is today's on the MCI-012 matrix", {
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

test_that("I5 to I13 the trx import is today's on the MCI-026 corpus and fixtures", {
    skip_if_not(workflowOwnerCapable())
    for (dir in c(module_ci_trx_corpus(), module_ci_trx_fixture("counts_whole_double"),
        module_ci_trx_fixture("counts_na"), module_ci_trx_fixture("counts_negative"),
        module_ci_trx_fixture("counts_fractional"))) {
        compareTrxImport(file.path(dir, "counts.tsv"), dir)
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

test_that("coerceOmicsNumericLikeColumns converts only numeric-like text and factor columns", {
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
    # The rule is metab's (metabImportValuesAreNumericLike,
    # coerceMetabImportWorkflowPayloadSamples).
    for (column in names(table)) {
        expect_identical(omicsImportValuesAreNumericLike(table[[column]]),
            metabImportValuesAreNumericLike(table[[column]]), info = column)
    }
    small <- table[c("id", "text", "factor")]
    payload <- coerceMetabImportWorkflowPayloadSamples(list(assayList = list(a = small),
        columnMapping = list(sample_columns = c("text", "factor"))))
    expect_identical(coerceOmicsNumericLikeColumns(small, c("text", "factor")), payload$assayList$a)
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
