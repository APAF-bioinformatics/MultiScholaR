module_ci_lipid_import_paths <- function() {
  root <- tempfile("module-ci-lipid-import-")
  paths <- list(
    base_dir = root,
    source_dir = file.path(root, "source"),
    raw_dir = file.path(root, "raw")
  )
  lapply(paths, dir.create, recursive = TRUE, showWarnings = FALSE)
  paths
}

module_ci_lipid_import_workflow <- function() {
  workflow <- new.env(parent = emptyenv())
  workflow$data_tbl <- NULL
  workflow$data_format <- NULL
  workflow$data_type <- NULL
  workflow$column_mapping <- NULL
  workflow$processing_log <- list(setup_import = NULL)
  workflow$tab_status <- list(setup_import = "incomplete", design_matrix = "disabled")
  workflow$workflow_types <- character()
  workflow$state_manager <- list(
    setWorkflowType = function(type) {
      workflow$workflow_types <- c(workflow$workflow_types, type)
      invisible(type)
    }
  )
  workflow
}

module_ci_lipid_sample_cols <- function() {
  c("WT_1", "WT_2", "KO_1", "KO_2")
}

module_ci_lipid_lipidsearch_assay <- function(include_lipid_name = TRUE,
                                              include_lipid_class = TRUE,
                                              duplicate_lipid = FALSE,
                                              extra_metadata = TRUE,
                                              zero_values = FALSE,
                                              non_numeric_samples = FALSE) {
  lipid_names <- if (isTRUE(duplicate_lipid)) {
    c("PC 34:1", "PC 34:1", "SM 34:1")
  } else {
    c("PC 34:1", "TG 52:3", "SM 34:1")
  }
  samples <- data.frame(
    WT_1 = c(10000, if (isTRUE(zero_values)) 0 else 30000, 8000),
    WT_2 = c(10200, 29800, 8200),
    KO_1 = c(21000, 30500, if (isTRUE(zero_values)) 0 else 16000),
    KO_2 = c(21400, 30200, 15800),
    check.names = FALSE
  )
  if (isTRUE(non_numeric_samples)) {
    samples$WT_1 <- c("ok", "bad", "values")
  }

  assay <- data.frame(
    LipidName = lipid_names,
    LipidClass = c("PC", "TG", "SM"),
    FattyAcid = c("16:0/18:1", "16:0/18:1/18:2", "18:1/16:0"),
    IonType = c("[M+H]+", "[M+NH4]+", "[M+H]+"),
    BaseRt = c(5.23, 12.34, 4.91),
    CalcMz = c(760.585, 874.785, 703.575),
    Grade = c("A", "A", "B"),
    samples,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (!isTRUE(include_lipid_name)) {
    assay$LipidName <- NULL
  }
  if (!isTRUE(include_lipid_class)) {
    assay$LipidClass <- NULL
  }
  if (isTRUE(extra_metadata)) {
    assay$BatchFlag <- c("keep", "keep", "review")
  }
  assay
}

module_ci_lipid_msdial_assay <- function(include_annotation = TRUE,
                                         duplicate_feature = FALSE,
                                         zero_missing = TRUE,
                                         mode = "Positive") {
  peak_ids <- if (isTRUE(duplicate_feature)) c("P001", "P001", "P003") else c("P001", "P002", "P003")
  assay <- data.frame(
    "Peak ID" = peak_ids,
    "RT (min)" = c(1.25, 2.5, 4.75),
    "Precursor m/z" = c(760.585, 874.785, 703.575),
    Adduct = c("[M+H]+", "[M+NH4]+", "[M-H]-"),
    Formula = c("C42H82NO8P", "C55H104NO6", "C39H79N2O6P"),
    Ontology = c("PC", "TG", "SM"),
    "Ion mode" = mode,
    "Total Score" = c(92.5, 85.2, 77.7),
    WT_1 = c(100, if (isTRUE(zero_missing)) 0 else 110, 120),
    WT_2 = c(101, 111, if (isTRUE(zero_missing)) NA_real_ else 121),
    KO_1 = c(200, 210, 220),
    KO_2 = c(201, 211, 221),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (isTRUE(include_annotation)) {
    assay <- data.frame(
      assay[1],
      Name = c("PC 34:1", "TG 52:3", "SM 34:1"),
      assay[-1],
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }
  assay
}

module_ci_lipid_custom_assay <- function(sample_prefix = c("WT", "sample"),
                                         numeric_as_character = FALSE) {
  sample_prefix <- match.arg(sample_prefix)
  sample_names <- if (identical(sample_prefix, "WT")) {
    module_ci_lipid_sample_cols()
  } else {
    c("Sample A", "Sample B", "QC-01", "QC-02")
  }
  samples <- data.frame(
    setNames(
      list(c(10, 11, 12), c(13, 14, 15), c(20, 21, 22), c(23, 24, 25)),
      sample_names
    ),
    check.names = FALSE
  )
  if (isTRUE(numeric_as_character)) {
    samples[] <- lapply(samples, as.character)
  }

  data.frame(
    lipid_id = c("L1", "L2", "L3"),
    annotation = c("PC", "TG", "SM"),
    samples,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

module_ci_lipid_write_table <- function(data, path) {
  utils::write.table(data, file = path, sep = "\t", row.names = FALSE, quote = FALSE)
  path
}

module_ci_lipid_write_malformed <- function(path) {
  writeLines(c("not\ta\tvalid", "\"unterminated"), path)
  path
}

module_ci_lipid_import_run <- function(assay,
                                       paths = module_ci_lipid_import_paths(),
                                       workflow = module_ci_lipid_import_workflow(),
                                       assay_name = "LCMS_Pos",
                                       assay_file = NULL,
                                       assay2_file = NULL,
                                       assay2_name = "",
                                       vendor_format = "lipidsearch",
                                       detected_format = vendor_format,
                                       lipid_id_col = "LipidName",
                                       annotation_col = "LipidClass",
                                       sample_cols = module_ci_lipid_sample_cols(),
                                       sanitize_names = FALSE,
                                       is_pattern = "^IS_") {
  if (is.null(assay_file)) {
    assay_file <- module_ci_lipid_write_table(
      assay,
      file.path(paths$raw_dir, paste0(gsub("[/\\\\]", "_", assay_name), ".tsv"))
    )
  }
  notifications <- list()
  removed_notifications <- character()
  error_messages <- character()

  result <- runLipidImportProcessing(
    workflowData = workflow,
    assay1Name = assay_name,
    assay1Data = assay,
    assay1File = assay_file,
    assay2File = assay2_file,
    assay2Name = assay2_name,
    vendorFormat = vendor_format,
    detectedFormat = detected_format,
    lipidIdCol = lipid_id_col,
    annotationCol = annotation_col,
    sampleColumns = sample_cols,
    isPattern = is_pattern,
    sanitizeNames = sanitize_names,
    experimentPaths = paths,
    notify = function(message, type = NULL, id = NULL, duration = NULL) {
      notifications[[length(notifications) + 1L]] <<- list(
        message = message,
        type = type,
        id = id,
        duration = duration
      )
      invisible(NULL)
    },
    removeNotify = function(id) {
      removed_notifications <<- c(removed_notifications, id)
      invisible(NULL)
    },
    logError = function(message) {
      error_messages <<- c(error_messages, message)
      invisible(NULL)
    }
  )

  list(
    result = result,
    workflow = workflow,
    paths = paths,
    notifications = notifications,
    removed_notifications = removed_notifications,
    error_messages = error_messages
  )
}

module_ci_lipid_preview <- function(path) {
  loadLipidImportAssayPreview(path)
}

module_ci_lipid_import_artifacts <- function(paths, assay_names) {
  c(
    setNames(file.path(paths$source_dir, paste0("data_cln_", assay_names, ".tab")), paste0("data_cln_", assay_names)),
    assay_manifest = file.path(paths$source_dir, "assay_manifest.txt"),
    column_mapping = file.path(paths$source_dir, "column_mapping.json"),
    import_summary = file.path(paths$source_dir, "lipidomics_import_summary.tsv")
  )
}

module_ci_lipid_import_state_digest <- function(workflow) {
  list(
    setup_import = workflow$tab_status$setup_import,
    design_matrix = workflow$tab_status$design_matrix,
    column_mapping = workflow$column_mapping,
    processing_log = workflow$processing_log$setup_import,
    workflow_types = workflow$workflow_types
  )
}
