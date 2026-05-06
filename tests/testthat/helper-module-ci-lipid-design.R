module_ci_lipid_design_paths <- function() {
  root <- tempfile("module-ci-lipid-design-")
  paths <- list(
    base_dir = root,
    source_dir = file.path(root, "source"),
    raw_dir = file.path(root, "raw")
  )
  lapply(paths, dir.create, recursive = TRUE, showWarnings = FALSE)
  paths
}

module_ci_lipid_design_samples <- function(kind = "balanced") {
  switch(
    kind,
    three_group = c("WT_1", "WT_2", "KO_1", "KO_2", "RES_1", "RES_2"),
    unbalanced = c("WT_1", "WT_2", "WT_3", "KO_1", "KO_2"),
    alternate = c("Sample-A", "Sample-B", "QC-01", "QC-02"),
    c("WT_1", "WT_2", "KO_1", "KO_2")
  )
}

module_ci_lipid_design_assay <- function(samples = module_ci_lipid_design_samples(),
                                         lipid_id_col = "LipidName",
                                         annotation_col = "LipidClass",
                                         lipid_prefix = "L",
                                         sample_order = samples) {
  sample_values <- setNames(
    lapply(seq_along(sample_order), function(idx) c(1000, 2000, 3000) + (idx * 100)),
    sample_order
  )
  data.frame(
    setNames(list(paste0(lipid_prefix, sprintf("%03d", 1:3))), lipid_id_col),
    setNames(list(c("PC", "TG", "SM")), annotation_col),
    IonMode = c("Positive", "Positive", "Negative"),
    sample_values,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

module_ci_lipid_design_assays <- function(layout = "lc_pair",
                                          samples = module_ci_lipid_design_samples(),
                                          reorder_secondary = FALSE,
                                          missing_sample = NULL,
                                          extra_sample = NULL,
                                          case_varied = FALSE) {
  assay_samples <- samples
  if (!is.null(missing_sample)) {
    assay_samples <- setdiff(assay_samples, missing_sample)
  }
  if (!is.null(extra_sample)) {
    assay_samples <- c(assay_samples, extra_sample)
  }
  if (isTRUE(case_varied)) {
    assay_samples[assay_samples == "WT_1"] <- "wt_1"
  }

  primary <- module_ci_lipid_design_assay(samples = samples, sample_order = assay_samples, lipid_prefix = "POS")
  secondary_order <- if (isTRUE(reorder_secondary)) rev(assay_samples) else assay_samples
  secondary <- module_ci_lipid_design_assay(samples = samples, sample_order = secondary_order, lipid_prefix = "NEG")

  switch(
    layout,
    gc = list(GCMS = primary),
    combined = list(LCMS_Pos = primary, GCMS = secondary),
    lc_pair = list(LCMS_Pos = primary, LCMS_Neg = secondary),
    list(LCMS_Pos = primary)
  )
}

module_ci_lipid_design_column_mapping <- function(samples = module_ci_lipid_design_samples(),
                                                  lipid_id_col = "LipidName",
                                                  annotation_col = "LipidClass") {
  list(
    lipid_id_col = lipid_id_col,
    annotation_col = annotation_col,
    sample_columns = samples,
    is_pattern = NA_character_
  )
}

module_ci_lipid_design_matrix <- function(kind = "two_group", samples = NULL) {
  if (is.null(samples)) {
    samples <- module_ci_lipid_design_samples(
      if (identical(kind, "three_group")) {
        "three_group"
      } else if (identical(kind, "unbalanced")) {
        "unbalanced"
      } else {
        "balanced"
      }
    )
  }
  group <- ifelse(grepl("^WT|Sample", samples), "WT", ifelse(grepl("^RES", samples), "RES", "KO"))
  batch <- rep(c("B1", "B2"), length.out = length(samples))

  design <- data.frame(
    Run = samples,
    factor1 = group,
    factor2 = if (identical(kind, "batch_aware")) batch else NA_character_,
    factor3 = NA_character_,
    group = if (identical(kind, "batch_aware")) paste(group, batch, sep = "_") else group,
    batch = batch,
    replicates = ave(seq_along(samples), group, FUN = seq_along),
    tech_reps = NA_integer_,
    stringsAsFactors = FALSE
  )

  if (identical(kind, "missing_group")) {
    design$group[[1]] <- NA_character_
  }
  if (identical(kind, "extra_metadata")) {
    design$site <- rep(c("left", "right"), length.out = nrow(design))
    design$instrument_batch <- rep(c("MS1", "MS2"), length.out = nrow(design))
  }
  design
}

module_ci_lipid_design_contrasts <- function(kind = "valid") {
  switch(
    kind,
    raw_terms = data.frame(contrast_string = "groupKO-groupWT", stringsAsFactors = FALSE),
    friendly = data.frame(
      contrasts = "groupKO-groupWT",
      friendly_names = "KO_vs_WT",
      full_format = "KO_vs_WT=groupKO-groupWT",
      stringsAsFactors = FALSE
    ),
    reversed = data.frame(contrasts = "groupWT-groupKO", stringsAsFactors = FALSE),
    duplicate = data.frame(contrasts = c("groupKO-groupWT", "groupKO-groupWT"), stringsAsFactors = FALSE),
    invalid_term = data.frame(contrasts = "groupMISSING-groupWT", stringsAsFactors = FALSE),
    empty = data.frame(contrasts = character(), stringsAsFactors = FALSE),
    no_contrasts = NULL,
    data.frame(contrasts = "groupKO-groupWT", stringsAsFactors = FALSE)
  )
}

module_ci_lipid_design_config <- function(formula = "~ 0 + group") {
  list(
    globalParameters = list(workflow_type = "lipidomics_standard"),
    deAnalysisParameters = list(formula_string = formula)
  )
}

module_ci_lipid_design_payload <- function(kind = "two_group",
                                           layout = "lc_pair",
                                           formula = "~ 0 + group",
                                           contrasts = module_ci_lipid_design_contrasts("friendly")) {
  design <- module_ci_lipid_design_matrix(kind)
  assays <- module_ci_lipid_design_assays(layout = layout, samples = design$Run)
  list(
    design_matrix = design,
    data_cln = assays,
    contrasts_tbl = contrasts,
    config_list = module_ci_lipid_design_config(formula),
    column_mapping = module_ci_lipid_design_column_mapping(design$Run)
  )
}

module_ci_lipid_design_write_import_pack <- function(root, payload) {
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  utils::write.table(payload$design_matrix, file.path(root, "design_matrix.tab"), sep = "\t", quote = FALSE, row.names = FALSE)
  for (assay_name in names(payload$data_cln)) {
    utils::write.table(payload$data_cln[[assay_name]], file.path(root, paste0("data_cln_", assay_name, ".tab")), sep = "\t", quote = FALSE, row.names = FALSE)
  }
  writeLines(names(payload$data_cln), file.path(root, "assay_manifest.txt"), useBytes = TRUE)
  if (!is.null(payload$contrasts_tbl) && nrow(payload$contrasts_tbl) > 0L) {
    writeLines(payload$contrasts_tbl$contrasts, file.path(root, "contrast_strings.tab"), useBytes = TRUE)
  }
  jsonlite::write_json(payload$column_mapping, file.path(root, "column_mapping.json"), auto_unbox = TRUE, pretty = TRUE, null = "null")
  writeLines("[globalParameters]\nworkflow_type = lipidomics_standard\n[deAnalysisParameters]\nformula_string = ~ 0 + group", file.path(root, "config.ini"), useBytes = TRUE)
  root
}

module_ci_lipid_design_workflow <- function(payload = NULL) {
  workflow <- new.env(parent = emptyenv())
  workflow$data_tbl <- if (is.null(payload)) NULL else payload$data_cln
  workflow$data_cln <- NULL
  workflow$design_matrix <- NULL
  workflow$contrasts_tbl <- NULL
  workflow$column_mapping <- if (is.null(payload)) NULL else payload$column_mapping
  workflow$config_list <- if (is.null(payload)) module_ci_lipid_design_config() else payload$config_list
  workflow$tab_status <- list(setup_import = "complete", design_matrix = "pending", quality_control = "disabled")
  workflow$state_manager <- NULL
  workflow
}

module_ci_lipid_design_assert_samples_once <- function(design, expected_samples) {
  testthat::expect_setequal(as.character(design$Run), expected_samples)
  testthat::expect_false(anyDuplicated(as.character(design$Run)) > 0L)
  invisible(TRUE)
}

module_ci_lipid_design_assert_preflight_passes <- function(payload, require_contrasts = TRUE) {
  result <- validateLipidDesignDaPreflight(
    designMatrix = payload$design_matrix,
    assayList = payload$data_cln,
    contrastsTbl = payload$contrasts_tbl,
    formulaString = payload$config_list$deAnalysisParameters$formula_string,
    columnMapping = payload$column_mapping,
    requireContrasts = require_contrasts
  )
  testthat::expect_true(result$valid, info = paste(result$errors, collapse = "; "))
  invisible(result)
}

module_ci_lipid_design_file_nonempty <- function(path) {
  testthat::expect_true(file.exists(path))
  testthat::expect_gt(file.info(path)$size, 0)
  invisible(TRUE)
}
