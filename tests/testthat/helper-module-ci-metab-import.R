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
  # KIT-072: the session's pending on-disk owner (the modules' only state manager) records the import's type; it is
  # closed when the test file ends.
  omicsTestOwner(workflow = workflow, omic_type = "metabolomics", bind = FALSE, env = testthat::teardown_env())
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
  notify <- function(message, ...) {
    notifications[[length(notifications) + 1L]] <<- list(message = message, ...)
    invisible(NULL)
  }
  quiet <- function(...) invisible(NULL)
  # KIT-102 W5: the import tab's processing is the shared feature-table panel's (the core import, then the commit and
  # the frame's completion), under the metabolomics configuration. Its result is `list(status, payload)` or
  # `list(status = "error", error)`.
  # KIT-075: the import binds the session's owner, which attaches to the reactive domain.
  config <- omicGuiConfig("metabolomics")
  result <- shiny::withReactiveDomain(omicsOwnerSession(), processOmicImportFeatureTable(
    workflow,
    request = list(assay1_name = assay_name, assay1_data = assay, assay1_file = assay_file, assay2_file = assay2_file,
      assay2_name = assay2_name, vendor_format = vendor_format, detected_format = detected_format,
      id_col = metabolite_col, annotation_col = annotation_col, sample_columns = sample_cols,
      sanitize_names = sanitize_names, is_pattern = is_pattern, experiment_paths = paths),
    config = config,
    complete = function(...) completeOmicImport(workflow, config, notify = FALSE),
    showNotification = notify,
    removeNotification = quiet,
    logError = quiet,
    logInfo = function(message) {
      logs <<- c(logs, message)
      invisible(NULL)
    }
  ))
  if (identical(result$status, "success")) result$artifactResult <- result$payload$artifacts
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
