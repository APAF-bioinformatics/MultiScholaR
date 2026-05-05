module_ci_metab_import_paths <- function() {
  root <- tempfile("module-ci-metab-import-")
  paths <- list(
    base_dir = root,
    source_dir = file.path(root, "source"),
    raw_dir = file.path(root, "raw")
  )
  lapply(paths, dir.create, recursive = TRUE, showWarnings = FALSE)
  paths
}

module_ci_metab_import_workflow <- function() {
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

module_ci_metab_custom_assay <- function(annotation = TRUE,
                                         numeric_as_character = FALSE,
                                         sample_prefix = c("WT", "KO"),
                                         extra_metadata = TRUE) {
  sample_prefix <- match.arg(sample_prefix)
  sample_names <- if (identical(sample_prefix, "WT")) {
    c("WT_1", "WT_2", "KO_1", "KO_2")
  } else {
    c("Sample-A", "Sample-B", "QC-01", "QC-02")
  }
  samples <- data.frame(
    setNames(
      list(c(1000, 1010, 1500), c(1005, 1008, 1495), c(1800, 1780, 1220), c(1810, 1775, 1215)),
      sample_names
    ),
    check.names = FALSE
  )
  if (isTRUE(numeric_as_character)) {
    samples[] <- lapply(samples, as.character)
  }

  assay <- data.frame(
    Feature.Name = c("M001", "M002", "M003"),
    samples,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (isTRUE(annotation)) {
    assay$Annotation <- c("Alanine", "Lactate", "Citrate")
  }
  if (isTRUE(extra_metadata)) {
    assay$Class <- c("Amino acid", "Organic acid", "Organic acid")
    assay$BatchFlag <- c("keep", "keep", "review")
  }
  assay
}

module_ci_metab_msdial_assay <- function(include_name = TRUE,
                                         duplicate_feature = FALSE,
                                         zero_missing = TRUE,
                                         mode = "Positive") {
  peak_ids <- if (isTRUE(duplicate_feature)) c("P001", "P001", "P003") else c("P001", "P002", "P003")
  assay <- data.frame(
    "Peak ID" = peak_ids,
    "RT (min)" = c(1.25, 2.5, 4.75),
    "Precursor m/z" = c(101.1, 202.2, 303.3),
    Adduct = c("[M+H]+", "[M+Na]+", "[M-H]-"),
    Formula = c("C3H7NO2", "C3H6O3", "C6H8O7"),
    Ontology = c("Amino acids", "Organic acids", "Organic acids"),
    "Ion mode" = mode,
    "Total Score" = c(92.5, 85.2, 77.7),
    WT_1 = c(100, if (isTRUE(zero_missing)) 0 else 110, 120),
    WT_2 = c(101, 111, if (isTRUE(zero_missing)) NA_real_ else 121),
    KO_1 = c(200, 210, 220),
    KO_2 = c(201, 211, 221),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (isTRUE(include_name)) {
    assay <- data.frame(
      assay[1],
      Name = c("Alanine", "Lactate", "Citrate"),
      assay[-1],
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }
  assay
}

module_ci_metab_write_table <- function(data, path) {
  utils::write.table(data, file = path, sep = "\t", row.names = FALSE, quote = FALSE)
  path
}

module_ci_metab_import_run <- function(assay,
                                       paths = module_ci_metab_import_paths(),
                                       workflow = module_ci_metab_import_workflow(),
                                       assay_name = "LCMS_Pos",
                                       assay_file = NULL,
                                       assay2_file = NULL,
                                       assay2_name = "",
                                       vendor_format = "custom",
                                       detected_format = "custom",
                                       metabolite_col = "Feature.Name",
                                       annotation_col = "Annotation",
                                       sample_cols = grep("^(WT|KO|Sample|QC)", names(assay), value = TRUE),
                                       sanitize_names = FALSE,
                                       is_pattern = "^IS_") {
  if (is.null(assay_file)) {
    assay_file <- module_ci_metab_write_table(
      assay,
      file.path(paths$raw_dir, paste0(gsub("[/\\\\]", "_", assay_name), ".tsv"))
    )
  }
  notifications <- list()
  logs <- character()
  result <- runMetabImportProcessing(
    assay1Data = assay,
    assay1Name = assay_name,
    assay1File = assay_file,
    assay2File = assay2_file,
    assay2Name = assay2_name,
    vendorFormat = vendor_format,
    detectedFormat = detected_format,
    sanitizeNames = sanitize_names,
    isPattern = is_pattern,
    getMetaboliteIdColFn = function() metabolite_col,
    getAnnotationColFn = function() annotation_col,
    getSampleColumnsFn = function() sample_cols,
    workflowData = workflow,
    experimentPaths = paths,
    reqFn = function(value) {
      if (is.null(value)) {
        stop("required value missing")
      }
      invisible(value)
    },
    showNotificationFn = function(message, ...) {
      notifications[[length(notifications) + 1L]] <<- list(message = message, ...)
      invisible(NULL)
    },
    logInfoFn = function(message) {
      logs <<- c(logs, message)
      invisible(NULL)
    },
    finalizeFeedbackFn = function(status, error = NULL) {
      invisible(list(status = status, error = if (is.null(error)) NULL else conditionMessage(error)))
    }
  )
  list(
    result = result,
    workflow = workflow,
    paths = paths,
    notifications = notifications,
    logs = logs
  )
}

module_ci_metab_import_artifacts <- function(paths, assay_names) {
  c(
    setNames(file.path(paths$source_dir, paste0("data_cln_", assay_names, ".tab")), paste0("data_cln_", assay_names)),
    assay_manifest = file.path(paths$source_dir, "assay_manifest.txt"),
    column_mapping = file.path(paths$source_dir, "column_mapping.json"),
    import_summary = file.path(paths$source_dir, "metabolomics_import_summary.tsv")
  )
}
