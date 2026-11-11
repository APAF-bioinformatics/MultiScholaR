# KIT-116 W2, quarantine row Q-30: the frozen values of test-workflow-omics-import.R, computed from the pre-KIT-102
# lipid, metab and trx import functions at their green commit d469cea3 (KIT-102 W2), where the comparison blocks
# passed. Run from a checkout of d469cea3 (git worktree add --detach <dir> d469cea3):
#
#   Rscript q30-import.R <checkout of d469cea3> <output dir>
#
# The script loads the package and the test helpers of that checkout, evaluates the test file's own helper
# definitions (compareLipidImport, compareMetabImport, compareTrxImport and the rest), and replaces only
# expectOmicsImportLikeToday by a recorder, so every case runs the block's own inputs through the old functions.
# Each value is recorded as the block compared it (the canonical forms below), and checked against the core of the
# same commit (the block's own comparison), so a digest stands for the identity the block asserted. Output:
# <output dir>/q30-frozen.R (the literals pasted into the test) and q30-old-values.rds (the old values themselves).

args <- commandArgs(trailingOnly = TRUE)
green <- normalizePath(args[[1]])
out_dir <- args[[2]]
stopifnot(identical(system2("git", c("-C", green, "rev-parse", "--short=8", "HEAD"), stdout = TRUE), "d469cea3"))
pkgload::load_all(green, quiet = TRUE)
withr::local_dir(green)
# The helpers defer their temporary owners to testthat's teardown environment, which test_file initialises.
testthat:::local_teardown_env()
testthat::source_test_helpers(file.path(green, "tests", "testthat"), env = globalenv())
library(testthat)
stopifnot(workflowOwnerCapable())

# The test file's top-level definitions (not its test_that blocks), as at the green commit.
exprs <- parse(file.path(green, "tests", "testthat", "test-workflow-omics-import.R"), keep.source = FALSE)
for (e in exprs) if (!(is.call(e) && identical(e[[1L]], as.name("test_that")))) eval(e, globalenv())

digestOf <- function(x) digest::digest(x, algo = "sha256")

# Each group as R code, one element per line; the code is checked to evaluate to the frozen value.
frozenCode <- function(frozen, prefix) {
    vapply(names(frozen), function(group) {
        values <- frozen[[group]]
        keys <- names(values)
        keys <- ifelse(make.names(keys) == keys, keys, paste0("`", keys, "`"))
        lines <- vapply(seq_along(values), function(i) {
            paste0("    ", keys[[i]], " = ", paste(deparse(values[[i]], width.cutoff = 500L), collapse = " "))
        }, character(1))
        code <- paste0(prefix, group, " <- list(\n", paste(lines, collapse = ",\n"), ")")
        stopifnot(identical(eval(parse(text = code)), values))
        code
    }, character(1))
}

# The canonical forms (the same functions are in the test file): a success as the block compared it, paths under the
# run's root made relative; a preview's fields under the core's names.
omicsImportRecord <- function(tables, sample_columns, data_format, column_mapping, log, validation, artifacts, root) {
    list(tables = tables, sample_columns = sample_columns, data_format = data_format, column_mapping = column_mapping,
        log = omicsImportLogWithoutTime(log, root), validation = validation[c("valid", "errors", "warnings")],
        artifacts = omicsImportArtifactBytes(artifacts, root))
}

frozen <- list()
old_values <- list()
faithful <- character()
record <- function(group, key, value, core_value) {
    # `[<-` with a list keeps a NULL value (a resolution can be NULL).
    if (is.null(frozen[[group]])) frozen[[group]] <<- list()
    frozen[[group]][key] <<- list(value)
    if (!identical(value, core_value)) faithful <<- c(faithful, paste(group, key))
}

# The recorder that replaces the comparison: the old side frozen, the core of the same commit checked against it.
expectOmicsImportLikeToday <- function(old, new, old_root, new_paths, info) {
    group <- current_group
    old_values[[group]][[info]] <<- old
    if (is.null(old$tables)) {
        core <- if (inherits(new, "error")) conditionMessage(new) else "<core succeeded>"
        if (length(list.files(new_paths$source_dir, recursive = TRUE)) > 0L) faithful <<- c(faithful, paste(info, "files"))
        record(group, info, list(message = old$message), list(message = core))
        return(invisible(NULL))
    }
    value <- omicsImportRecord(old$tables, old$sample_columns, old$data_format, old$column_mapping, old$log,
        old$validation, old$artifacts, old_root)
    core <- if (inherits(new, "error")) NULL else omicsImportRecord(new$assays, new$sample_columns, new$data_format,
        new$column_mapping, new$processing_log, new$validation, new$artifacts, new_paths$base_dir)
    record(group, info, list(digest = digestOf(value)), list(digest = digestOf(core)))
    if (!identical(value, core)) faithful <<- c(faithful, paste(group, info, "value"))
}

# I1: the lipid and metab default header readers.
lipid_default <- eval(formals(loadLipidImportAssayPreview)$readHeaders)
metab_default <- eval(formals(prepareMetabImportAssaySelectionState)$readHeadersFn)
paths <- module_ci_lipid_import_paths()
files <- c(
    tsv = module_ci_lipid_write_table(module_ci_lipid_msdial_assay(), file.path(paths$raw_dir, "a.tsv")),
    malformed = module_ci_lipid_write_malformed(file.path(paths$raw_dir, "m.tsv")),
    missing = file.path(paths$raw_dir, "absent.tsv"))
csv <- file.path(paths$raw_dir, "q.csv")
utils::write.csv(module_ci_lipid_custom_assay(), csv, row.names = FALSE)
files <- c(files, csv = csv)
for (name in names(files)) {
    lipid <- suppressWarnings(lipid_default(files[[name]]))
    metab <- suppressWarnings(metab_default(files[[name]]))
    core <- suppressWarnings(readOmicsImportHeaders(files[[name]]))
    record("i1_lipid", name, lipid, core)
    record("i1_metab", name, metab, core)
}

# I2: the lipid and metab previews, as the fields the block compared, paths relative to the file's root.
previewRecord <- function(headers, detected_format, format_confidence, import_result, id, annotation, is_pattern,
                          assay_data, root) {
    omicsImportRelative(list(headers = headers, detected_format = detected_format,
        format_confidence = format_confidence, import_result = import_result, assay_data = assay_data, id = id,
        annotation = annotation, is_pattern = is_pattern), root)
}
lipid_files <- c(
    lipidsearch = module_ci_lipid_write_table(module_ci_lipid_lipidsearch_assay(), file.path(paths$raw_dir, "ls.txt")),
    msdial = module_ci_lipid_write_table(module_ci_lipid_msdial_assay(),
        file.path(paths$raw_dir, "msdial_positive.tsv")),
    no_annotation = module_ci_lipid_write_table(module_ci_lipid_msdial_assay(include_annotation = FALSE),
        file.path(paths$raw_dir, "msdial_no_annotation.tsv")),
    custom = module_ci_lipid_write_table(module_ci_lipid_custom_assay(), file.path(paths$raw_dir, "custom.tsv")),
    malformed = module_ci_lipid_write_malformed(file.path(paths$raw_dir, "malformed.tsv")))
lipid_config <- omicGuiConfig("lipidomics")$import
for (name in names(lipid_files)) {
    old <- suppressWarnings(loadLipidImportAssayPreview(lipid_files[[name]], logInfo = function(...) NULL))
    new <- suppressWarnings(previewOmicsImport(lipid_files[[name]], lipid_config))
    old_values$i2_lipid[[name]] <- old
    value <- previewRecord(old$headers, old$detectedFormat, old$formatConfidence, old$importResult,
        old$updates$lipidId, old$updates$annotation, old$updates$isPattern, old$assayData, paths$base_dir)
    core <- previewRecord(new$headers, new$detected_format, new$format_confidence, new$import_result,
        new$updates$id, new$updates$annotation, new$updates$is_pattern, new$assay_data, paths$base_dir)
    record("i2_lipid", name, digestOf(value), digestOf(core))
}
metab_config <- omicGuiConfig("metabolomics")$import
metab_paths <- module_ci_metab_import_paths()
metab_files <- c(
    msdial = module_ci_metab_write_table(module_ci_metab_msdial_assay(),
        file.path(metab_paths$raw_dir, "msdial_positive.tsv")),
    custom = module_ci_metab_write_table(module_ci_metab_custom_assay(), file.path(metab_paths$raw_dir, "custom.tsv")))
for (name in names(metab_files)) {
    old <- suppressWarnings(prepareMetabImportAssaySelectionState(metab_files[[name]]))
    new <- suppressWarnings(previewOmicsImport(metab_files[[name]], metab_config))
    old_values$i2_metab[[name]] <- old
    # The block compared these fields only (not the assay data).
    value <- previewRecord(old$headers, old$formatInfo$format, old$formatInfo$confidence, old$importResult,
        list(choices = old$metaboliteIdChoices, selected = old$selectedMetaboliteId),
        list(choices = old$annotationChoices, selected = old$selectedAnnotation), old$isPattern, NULL,
        metab_paths$base_dir)
    core <- previewRecord(new$headers, new$detected_format, new$format_confidence, new$import_result,
        new$updates$id, new$updates$annotation, new$updates$is_pattern, NULL, metab_paths$base_dir)
    record("i2_metab", name, digestOf(value), digestOf(core))
}

# I2 trx: trxImportPreview on the corpus and the negative fixture (the directory's name is the key).
for (dir in c(module_ci_trx_corpus(), module_ci_trx_fixture("counts_negative"))) {
    old <- trxImportPreview(file.path(dir, "counts.tsv"))
    new <- previewTrxCountsImport(file.path(dir, "counts.tsv"))
    old_values$i2_trx[[basename(dir)]] <- old
    value <- previewRecord(old$headers, old$detectedFormat, old$formatConfidence, old$importResult,
        old$updates$lipidId, old$updates$annotation, NULL, old$assayData, dir)
    core <- previewRecord(new$headers, new$detected_format, new$format_confidence, new$import_result,
        new$updates$id, new$updates$annotation, NULL, new$assay_data, dir)
    record("i2_trx", basename(dir), digestOf(value), digestOf(core))
}

# I3, I4: the sample and column resolution (small values, frozen as literals).
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
        record("i3_lipid", paste(i, exclude),
            resolveLipidImportSampleColumns(custom, case$result, case$vendor, case$pattern %||% "",
                excludeNormalized = exclude),
            resolveOmicsImportSampleColumns(custom, case$vendor, case$pattern, case$result,
                if (exclude) OMICS_IMPORT_NORMALISED_COLUMNS else NULL))
    }
    record("i3_metab", as.character(i), resolveMetabImportSampleColumns(custom, case$vendor, case$pattern, case$result),
        resolveOmicsImportSampleColumns(custom, case$vendor, case$pattern, case$result))
}
headers <- names(custom)
for (vendor in c("custom", "msdial")) {
    for (custom_name in list("LIPID_ID", "absent", "", NULL)) {
        key <- paste(vendor, custom_name %||% "NULL")
        new <- resolveOmicsImportColumn(headers, vendor, "selected", custom_name)
        record("i4_lipid", key, resolveLipidImportEffectiveColumn(custom, vendor, "selected", custom_name), new)
        metab_today <- if (identical(vendor, "custom")) resolveMetabImportColumnName(headers, custom_name) else "selected"
        record("i4_metab", key, metab_today, new)
    }
}

# I5 to I13, lipid: the block's own cases through compareLipidImport (the recorder above).
current_group <- "i5_lipid"
paths <- module_ci_lipid_import_paths()
for (case in list(
    list("ls_pos", module_ci_lipid_lipidsearch_assay()),
    list("ls_zero", module_ci_lipid_lipidsearch_assay(zero_values = TRUE), assay_name = "LCMS_Neg"),
    list("ls_no_class", module_ci_lipid_lipidsearch_assay(include_lipid_class = FALSE), annotation_col = ""),
    list("ls_duplicate_ids", module_ci_lipid_lipidsearch_assay(duplicate_lipid = TRUE)),
    list("missing_id", module_ci_lipid_lipidsearch_assay(include_lipid_name = FALSE)),
    list("missing_id_column", module_ci_lipid_lipidsearch_assay(), lipid_id_col = "MissingLipid"),
    list("no_samples", module_ci_lipid_lipidsearch_assay(), sample_cols = character()),
    list("non_numeric", module_ci_lipid_lipidsearch_assay(non_numeric_samples = TRUE)),
    list("mismatched", module_ci_lipid_lipidsearch_assay(), sample_cols = c("Missing_1", "Missing_2")),
    list("no_numeric", data.frame(LipidName = c("L1", "L2"), WT_1 = c("low", "high")), sample_cols = "WT_1"),
    list("msdial", module_ci_lipid_msdial_assay(duplicate_feature = TRUE), vendor_format = "msdial",
        lipid_id_col = "Peak ID", annotation_col = "Name"),
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
samples <- c("Sample A", "Sample B", "QC-01", "QC-02")
custom <- module_ci_lipid_custom_assay(sample_prefix = "sample")
assays <- list(a = custom, b = custom[c("lipid_id", "Sample A", "QC-02")])
for (sanitize in c(TRUE, FALSE)) {
    old <- sanitizeLipidImportSampleNames(assays, samples, sanitize, logInfo = function(...) NULL,
        notify = function(...) NULL)
    new <- sanitizeOmicsImportSampleNames(assays, samples, sanitize)
    record("i5_sanitize", as.character(sanitize),
        list(assays = digestOf(old$assayList), sample_columns = old$sampleColumns),
        list(assays = digestOf(new$assays), sample_columns = new$sample_columns))
}
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
msdial_neg <- module_ci_lipid_write_table(module_ci_lipid_msdial_assay(), file.path(paths$raw_dir, "msdial_neg.tsv"))
compareLipidImport(primary, "pair_msdial_reader", vendor_format = "msdial", detected_format = "msdial",
    assay2_file = msdial_neg, assay2_name = "LCMS_Neg")

# I5 to I13, metab.
current_group <- "i5_metab"
paths <- module_ci_metab_import_paths()
valid <- module_ci_metab_custom_assay()
for (case in list(
    list("annotation_present", module_ci_metab_custom_assay(numeric_as_character = FALSE)),
    list("annotation_absent", module_ci_metab_custom_assay(annotation = FALSE), annotation_col = ""),
    list("text_numeric", module_ci_metab_custom_assay(numeric_as_character = TRUE, sample_prefix = "KO"),
        sample_cols = c("Sample-A", "Sample-B", "QC-01", "QC-02")),
    list("msdial", module_ci_metab_msdial_assay(duplicate_feature = TRUE), vendor_format = "msdial",
        metabolite_col = "Peak ID", annotation_col = "Name", sample_cols = c("WT_1", "WT_2", "KO_1", "KO_2")),
    list("missing_id_column", valid, metabolite_col = "Missing.ID", sample_cols = c("WT_1", "WT_2")),
    list("empty_id", valid, metabolite_col = "", sample_cols = c("WT_1", "WT_2")),
    list("no_samples", valid, sample_cols = character()),
    list("non_numeric", within(valid, WT_1 <- c("ok", "bad", "values")), sample_cols = c("WT_1", "WT_2")),
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
secondary <- module_ci_metab_write_table(module_ci_metab_custom_assay(), file.path(paths$raw_dir, "lcms_neg.tsv"))
compareMetabImport(valid, "pair", assay2_file = secondary, assay2_name = "LCMS_Neg")
compareMetabImport(valid, "pair_unusual", assay_name = "GCMS-2026_A", assay2_file = secondary,
    assay2_name = "LCMS_Pos_2026")
compareMetabImport(valid, "pair_duplicate", assay2_file = secondary, assay2_name = "LCMS_Pos")
compareMetabImport(valid, "pair_no_name", assay_name = "GCMS", assay2_file = secondary, assay2_name = "")
compareMetabImport(valid, "pair_path", assay2_file = secondary, assay2_name = "LCMS\nNeg")

# I5 to I13, trx (the directory's name is the key).
current_group <- "i5_trx"
for (dir in c(module_ci_trx_corpus(), module_ci_trx_fixture("counts_whole_double"), module_ci_trx_fixture("counts_na"),
    module_ci_trx_fixture("counts_negative"), module_ci_trx_fixture("counts_fractional"))) {
    compareTrxImport(file.path(dir, "counts.tsv"), basename(dir))
}

# The metab cross-check of the shared conversion.
table <- data.frame(id = c("a", "b", "c"), text = c(" 1.5", "", NA), factor = factor(c("2", "3", "4")),
    word = c("1", "x", "2"), number = c(1, 2, 3), flag = c(NA, NA, NA), unlisted = c("7", "8", "9"),
    infinite = c("1", "Inf", "2"), stringsAsFactors = FALSE)
numeric_like <- vapply(table, metabImportValuesAreNumericLike, logical(1))
record("coerce", "numeric_like", numeric_like, vapply(table, omicsImportValuesAreNumericLike, logical(1)))
small <- table[c("id", "text", "factor")]
payload <- coerceMetabImportWorkflowPayloadSamples(list(assayList = list(a = small),
    columnMapping = list(sample_columns = c("text", "factor"))))
record("coerce", "payload", payload$assayList$a, coerceOmicsNumericLikeColumns(small, c("text", "factor")))

if (length(faithful)) stop("the old value is not the green core's for: ", paste(faithful, collapse = ", "))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(old_values, file.path(out_dir, "q30-old-values.rds"))
writeLines(frozenCode(frozen, "omicsImportFrozen_"), file.path(out_dir, "q30-frozen.R"))
cat("frozen groups:", paste(names(frozen), lengths(frozen), sep = "=", collapse = " "), "\n")
