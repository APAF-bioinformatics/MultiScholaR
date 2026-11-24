# KIT-198 gate G0: does a numeric column that is not a sample reach the lipid or metabolite intensity filter on the
# application's own import path? (ticket section 9 "G0", acceptance criterion A0; owner decision D8, methods ledger
# VAL-118 (1)).
#
# Usage, from the package root of the gate commit on a run-tests runner:
#   Rscript docs/qa/workflow-execution/kit/kit-198/g0-numeric-columns.R [<out dir>]
#   <out dir>  default docs/qa/workflow-execution/kit/kit-198/g0
#
# Under pkgload::load_all(".", export_all = TRUE, helpers = TRUE), each input is taken through the functions the import
# panel and the design builder call (R/mod_omic_import_feature_table.R, R/mod_omic_design_builder.R):
#   the configuration's preview (`previewOmicsImport`) for the detected format and the proposed id, annotation and
#   internal-standard columns (the values the panel's update* calls send back); the panel's column resolvers
#   (`resolveOmicsImportColumn`, `resolveOmicsImportSampleColumns`, the lipid panel's "Exclude Vendor-Normalized
#   Columns" default TRUE); `importOmicsTables(request, omicGuiConfig(<omic>))`; then
#   `buildOmicsDesignSaveDataList()` with the input's design; then the configuration's `design$build_object()`, whose
#   assay tables are what the intensity filter receives.
# Inputs:
#   the four GUI-entry lanes with their manifest import settings and the lane driver's inputs
#   (tests/testthat/helper-gui-entry.R guiEntryOmicsSteps: custom lanes id Feature.Name, no annotation, samples
#   "^(WT|KO)_", internal standards "^ITSD|_d[0-9]+$"; sanitising off) and design (design.tsv);
#   two MS-DIAL tables with integer ids: tests/testdata/e2e/metab_lc/seed_lcms_pos.csv through the metabolomics
#   MS-DIAL import, and module_ci_lipid_msdial_assay() with `Peak ID` = seq_len(nrow(.)) through the lipid MS-DIAL
#   import (panel defaults otherwise; every sample column the panel resolves is assigned except the id and annotation
#   columns, group = the name before the first "_"). The id and annotation columns are never counted as samples.
# Written:
#   g0-numeric-columns.tsv: input, assay, id column and its class, annotation column and its class, the design's
#     samples, the columns the filter counts (`names(assay_table)[sapply(assay_table, is.numeric)]`, as
#     R/func_lipid_qc.R:47 and R/func_metab_qc_filtering_helpers.R:20), the numeric columns that are not samples of
#     the design, and whether the design object's assay tables count the same columns;
#   g0-threshold-effect.tsv (reported only; no decision rests on it): per input and assay, the threshold at
#     percentiles 1 and 10 over the counted columns (as the S4 kernels, R/func_lipid_s4_normalization_methods.R:931-941,
#     R/func_metab_s4_qc_methods.R:84-94) and over the design's samples only, and, at proportion 0.5 and the 10th
#     percentile, the features removed under each (the helper on the table as it is; the helper on the table without
#     the non-sample numeric columns);
#   g0-decision.txt: "clean" when no input has a non-sample numeric column, else "found" with the inputs.
# Synthetic and repository fixtures only. Exit status non-zero when an input fails to run.

args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1L) args[[1L]] else "docs/qa/workflow-execution/kit/kit-198/g0"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out)
commit <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
cat("KIT-198 G0 at", commit, "out", out, "\n")

suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE, export_all = TRUE, helpers = TRUE))
ns <- asNamespace("MultiScholaR")
helper_fn <- function(name) get(name, envir = ns, inherits = TRUE)

e2e <- normalizePath(file.path("tests", "testdata", "e2e"), mustWork = TRUE)
manifest <- jsonlite::fromJSON(file.path(e2e, "manifest.json"), simplifyVector = FALSE)
lane_of <- function(id) {
    for (lane in manifest$lanes) if (identical(lane$lane_id, id)) return(lane)
    stop("lane not in the manifest: ", id)
}

# One input through the panel's request, the import, the design save and the design object.
run_input <- function(input_id, omic, vendor_format, assay1_file, assay1_name, assay2_file = NULL,
                      assay2_name = NULL, custom = NULL, sanitize = TRUE, design = NULL) {
    config <- omicGuiConfig(omic)
    import <- config$import
    preview <- config$import$preview(assay1_file, import = import)
    headers <- names(preview$assay_data)
    if (identical(vendor_format, "custom")) {
        id_col <- resolveOmicsImportColumn(headers, vendor_format, NULL, custom$id)
        annotation_col <- resolveOmicsImportColumn(headers, vendor_format, NULL, custom$annotation)
        sample_pattern <- custom$samples
        is_pattern <- custom$is_pattern
    } else {
        id_col <- preview$updates$id$selected
        annotation_col <- preview$updates$annotation$selected
        sample_pattern <- ""
        is_pattern <- preview$updates$is_pattern %||% ""
    }
    exclude <- if (identical(omic, "lipidomics")) OMICS_IMPORT_NORMALISED_COLUMNS
    sample_columns <- resolveOmicsImportSampleColumns(preview$assay_data, vendor_format, sample_pattern,
        preview$import_result, exclude)
    source_dir <- file.path(tempfile("kit198-g0-"), "source")
    dir.create(source_dir, recursive = TRUE)
    request <- list(assay1_name = assay1_name, assay1_data = preview$assay_data, assay1_file = assay1_file,
        assay2_file = assay2_file, assay2_name = assay2_name, vendor_format = vendor_format,
        detected_format = preview$detected_format, id_col = id_col, annotation_col = annotation_col,
        sample_columns = sample_columns, sanitize_names = sanitize, is_pattern = is_pattern,
        experiment_paths = list(source_dir = source_dir))
    imported <- importOmicsTables(request, config)
    if (is.null(design)) {
        id_and_annotation <- c(imported$column_mapping[[import$id_mapping_key]], imported$column_mapping$annotation_col)
        samples <- setdiff(imported$sample_columns, id_and_annotation)
        design <- data.frame(sample = samples, group = sub("_.*$", "", samples), stringsAsFactors = FALSE)
    }
    id_key <- import$id_mapping_key
    data_cln <- buildOmicsDesignSaveDataList(imported$assays, design$sample, imported$column_mapping, id_key)
    design_matrix <- data.frame(Run = design$sample, group = design$group,
        tech_rep_group = paste(design$group, ave(design$sample, design$group, FUN = seq_along), sep = "_"),
        stringsAsFactors = FALSE)
    object <- config$design$build_object(list(data_cln = data_cln, design_matrix = design_matrix,
        config_list = list()), imported$column_mapping)
    object_assays <- methods::slot(object, if (identical(omic, "lipidomics")) "lipid_data" else "metabolite_data")
    list(omic = omic, id_col = imported$column_mapping[[id_key]], annotation_col = imported$column_mapping$annotation_col,
        samples = design$sample, panel_samples = imported$sample_columns, data_cln = data_cln, object_assays = object_assays,
        detected_format = preview$detected_format)
}

numeric_columns <- function(table) names(table)[vapply(table, is.numeric, logical(1))]
join <- function(x) if (length(x)) paste(x, collapse = ";") else "none"
column_class <- function(table, column) {
    if (is.null(column) || !nzchar(column) || !column %in% names(table)) return("absent")
    paste(class(table[[column]]), collapse = "/")
}
threshold <- function(table, columns, percentile) {
    values <- unlist(table[, columns, drop = FALSE], use.names = FALSE)
    unname(stats::quantile(values, na.rm = TRUE, probs = percentile / 100)[1])
}
removed <- function(omic, table, threshold, proportion, id_col) {
    kept <- if (identical(omic, "lipidomics")) {
        lipidIntensityFilteringHelper(assay_table = table, min_lipid_intensity_threshold = threshold,
            lipids_proportion_of_samples_below_cutoff = proportion, lipid_id_column = id_col)
    } else {
        metaboliteIntensityFilteringHelper(assay_table = table, min_metabolite_intensity_threshold = threshold,
            metabolites_proportion_of_samples_below_cutoff = proportion, metabolite_id_column = id_col)
    }
    nrow(table) - nrow(kept)
}

custom_inputs <- list(id = "Feature.Name", annotation = "", samples = "^(WT|KO)_", is_pattern = "^ITSD|_d[0-9]+$")
lane_input <- function(lane_id) {
    lane <- lane_of(lane_id)
    dir <- file.path(e2e, lane$fixture_dir)
    assays <- unlist(lane$assays)
    design <- utils::read.delim(file.path(dir, "design.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
    second <- !is.null(lane$assay2_file) && length(assays) >= 2L
    function() run_input(lane_id, lane$omic_type, lane$import_tool, file.path(dir, lane$seed_file), assays[[1L]],
        assay2_file = if (second) file.path(dir, lane$assay2_file), assay2_name = if (second) assays[[2L]],
        custom = if (identical(lane$import_tool, "custom")) custom_inputs, sanitize = FALSE, design = design)
}
lipid_msdial <- function() {
    assay <- helper_fn("module_ci_lipid_msdial_assay")()
    assay[["Peak ID"]] <- seq_len(nrow(assay))
    path <- file.path(tempfile("kit198-g0-lipid-"), "msdial_lipid_pos.csv")
    dir.create(dirname(path), recursive = TRUE)
    utils::write.csv(assay, path, row.names = FALSE)
    run_input("msdial_lipid_module_ci", "lipidomics", "msdial", path, "LCMS_Pos")
}
metab_msdial <- function() {
    run_input("msdial_metab_seed_lcms_pos", "metabolomics", "msdial", file.path(e2e, "metab_lc", "seed_lcms_pos.csv"),
        "LCMS_Pos")
}
inputs <- list(lipid_canonical = lane_input("lipid_canonical"), metab_lc = lane_input("metab_lc"),
    metab_gc = lane_input("metab_gc"), metab_combined = lane_input("metab_combined"),
    msdial_lipid_module_ci = lipid_msdial, msdial_metab_seed_lcms_pos = metab_msdial)

columns_rows <- list()
effect_rows <- list()
failed <- character()
for (input_id in names(inputs)) {
    cat("input", input_id, "\n")
    result <- tryCatch(inputs[[input_id]](), error = function(e) {
        cat("  FAILED:", conditionMessage(e), "\n")
        failed <<- c(failed, input_id)
        NULL
    })
    if (is.null(result)) next
    for (assay_name in names(result$data_cln)) {
        table <- result$data_cln[[assay_name]]
        counted <- numeric_columns(table)
        # The id and annotation columns are never samples, even when the panel proposes one as a sample column.
        samples <- setdiff(intersect(result$samples, names(table)), c(result$id_col, result$annotation_col))
        non_sample <- setdiff(counted, samples)
        object_counted <- numeric_columns(result$object_assays[[assay_name]])
        columns_rows[[length(columns_rows) + 1L]] <- data.frame(input = input_id, omic = result$omic,
            detected_format = result$detected_format %||% "", assay = assay_name, id_column = result$id_col,
            id_class = column_class(table, result$id_col), annotation_column = result$annotation_col %||% "",
            annotation_class = column_class(table, result$annotation_col), n_features = nrow(table),
            panel_sample_columns = join(result$panel_samples), design_samples = join(samples), counted_columns = join(counted), n_counted = length(counted),
            n_samples = length(samples), non_sample_numeric = join(non_sample),
            object_counts_same = identical(object_counted, counted), stringsAsFactors = FALSE)
        sample_table <- table[, setdiff(names(table), non_sample), drop = FALSE]
        t_counted <- vapply(c(1, 10), \(p) threshold(table, counted, p), numeric(1))
        t_samples <- vapply(c(1, 10), \(p) threshold(table, samples, p), numeric(1))
        effect_rows[[length(effect_rows) + 1L]] <- data.frame(input = input_id, assay = assay_name,
            counted_columns = join(counted), non_sample_numeric = join(non_sample), n_features = nrow(table),
            threshold_p1_counted = t_counted[[1]], threshold_p1_samples_only = t_samples[[1]],
            threshold_p10_counted = t_counted[[2]], threshold_p10_samples_only = t_samples[[2]],
            removed_p10_prop0.5_counted = removed(result$omic, table, t_counted[[2]], 0.5, result$id_col),
            removed_p10_prop0.5_samples_only = removed(result$omic, sample_table, t_samples[[2]], 0.5, result$id_col),
            stringsAsFactors = FALSE)
    }
}

columns <- do.call(rbind, columns_rows)
effect <- do.call(rbind, effect_rows)
utils::write.table(columns, file.path(out, "g0-numeric-columns.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
utils::write.table(effect, file.path(out, "g0-threshold-effect.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
found <- unique(columns$input[columns$non_sample_numeric != "none"])
decision <- c(paste("commit", commit),
    paste("inputs run:", paste(unique(columns$input), collapse = ", ")),
    paste("inputs failed:", if (length(failed)) paste(failed, collapse = ", ") else "none"),
    if (length(found)) paste("found: a non-sample numeric column reaches the filter on", paste(found, collapse = ", "))
    else "clean: no input has a non-sample numeric column")
writeLines(decision, file.path(out, "g0-decision.txt"))
cat(decision, sep = "\n")
print(columns)
print(effect)
if (length(failed)) quit(status = 1L)
