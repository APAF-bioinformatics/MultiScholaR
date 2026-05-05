module_ci_prot_write_tsv <- function(data, path) {
  utils::write.table(
    data,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )
  path
}

module_ci_prot_temp_path <- function(stem, ext = ".tsv") {
  tempfile(pattern = paste0("module-ci-prot-", stem, "-"), fileext = ext)
}

module_ci_prot_make_fasta <- function(path = module_ci_prot_temp_path("fasta", ".fasta")) {
  writeLines(
    c(
      ">sp|P001|PROT1 OS=Homo sapiens OX=9606",
      "MPEPTIDESEQ",
      ">sp|P002|PROT2 OS=Homo sapiens OX=9606",
      "MPEPTIDESEQAA"
    ),
    path,
    useBytes = TRUE
  )
  path
}

module_ci_prot_make_diann <- function(path = module_ci_prot_temp_path("diann"), variant = "canonical") {
  data <- data.frame(
    Protein.Group = c("P001", "P001", "P002", "P003"),
    Protein.Ids = c("P001", "P001", "P002", "P003"),
    Protein.Names = c("Protein 1", "Protein 1", "Protein 2", "Protein 3"),
    Precursor.Id = c("P001_PEP1_2", "P001_PEP2_2", "P002_PEP1_2", "P003_PEP1_3"),
    Modified.Sequence = c("PEP1", "PEP2", "PEP3", "PEP4"),
    Stripped.Sequence = c("PEP1", "PEP2", "PEP3", "PEP4"),
    Precursor.Charge = c(2L, 2L, 2L, 3L),
    Q.Value = c(0.001, 0.002, 0.003, 0.004),
    PG.Q.Value = c(0.001, 0.001, 0.002, 0.003),
    Run = c("WT_1", "WT_2", "KO_1", "KO_2"),
    Precursor.Quantity = c(1000, 1200, 2500, 2600),
    Precursor.Normalised = c(10, 11, 20, 21),
    Extra.Metadata = c("keep_a", "keep_b", "keep_c", "keep_d"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (identical(variant, "duplicate_proteins")) {
    data$Protein.Group <- c("P001", "P001", "P001", "P002")
    data$Protein.Ids <- data$Protein.Group
  }
  if (identical(variant, "zero_na")) {
    data$Precursor.Quantity <- c(0, NA, 2500, 2600)
    data$Precursor.Normalised <- c(0, NA, 20, 21)
  }
  if (identical(variant, "bad_run_names")) {
    data$Run <- c("123 Sample-A", "WT sample 2", "KO/sample 1", "KO sample 2")
  }
  if (identical(variant, "missing_required")) {
    data$Run <- NULL
  }

  module_ci_prot_write_tsv(data, path)
}

module_ci_prot_make_tmt <- function(path = module_ci_prot_temp_path("tmt"), variant = "canonical") {
  if (identical(variant, "missing_accession")) {
    data <- data.frame(
      Annotated.Sequence = c("PEP1", "PEP2"),
      `Abundance: F1: 126, WT_1` = c(100, 200),
      check.names = FALSE
    )
    return(module_ci_prot_write_tsv(data, path))
  }

  data <- data.frame(
    Annotated.Sequence = c("PEP1", "PEP2", "PEP_SHARED"),
    Master.Protein.Accessions = c("P001", "P002", "P001;P002"),
    `Abundance: F1: 126, WT_1` = c(100, 200, 50),
    `Abundance: F2: 127N, WT_2` = c(110, 210, 55),
    `Abundance: F3: 128C, KO_1` = c(300, 400, 75),
    `Abundance: F4: 129N, KO_2` = c(310, 410, 80),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (identical(variant, "simplified")) {
    names(data) <- c(
      "Annotated.Sequence",
      "Master.Protein.Accessions",
      "Abundance.126",
      "Abundance.127N",
      "Abundance.128C",
      "Abundance.129N"
    )
  }
  if (identical(variant, "collision")) {
    names(data)[3:4] <- c("Abundance: F1: 126, WT_1", "Abundance: F2: 126, WT_1")
  }
  if (identical(variant, "all_na_channel")) {
    data[["Abundance: F4: 129N, KO_2"]] <- NA_real_
  }

  module_ci_prot_write_tsv(data, path)
}

module_ci_prot_make_maxquant <- function(path = module_ci_prot_temp_path("maxquant", ".txt"), variant = "canonical") {
  data <- data.frame(
    Protein.IDs = c("P001", "P002;P003", "CON__P004", "REV__P005"),
    Gene.names = c("G1", "", "GC", "GR"),
    Protein.names = c("Protein 1", "Protein group", "Contaminant", "Reverse"),
    LFQ.intensity.WT_1 = c(0, 2000, 3000, 4000),
    LFQ.intensity.WT_2 = c(100, 2100, 3100, 4100),
    LFQ.intensity.KO_1 = c(500, 2500, 3500, 4500),
    LFQ.intensity.KO_2 = c(600, 2600, 3600, 4600),
    Potential.contaminant = c("", "", "+", ""),
    Reverse = c("", "", "", "+"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (identical(variant, "duplicate_accessions")) {
    data$Protein.IDs[2] <- "P001"
  }
  if (identical(variant, "missing_intensity")) {
    data <- data[, c("Protein.IDs", "Gene.names"), drop = FALSE]
  }
  if (identical(variant, "missing_protein")) {
    data$Protein.IDs <- NULL
  }

  module_ci_prot_write_tsv(data, path)
}

module_ci_prot_make_fragpipe <- function(path = module_ci_prot_temp_path("fragpipe"), variant = "canonical") {
  data <- data.frame(
    Protein = c("P001", "P002"),
    `Protein ID` = c("sp|P001|PROT1", "sp|P002|PROT2"),
    Gene = c("G1", "G2"),
    `WT A MaxLFQ Intensity` = c(1000, 2000),
    `WT-B MaxLFQ Intensity` = c(1100, 2100),
    `KO/sample 1 MaxLFQ Intensity` = c(3000, 4000),
    `KO sample 2 MaxLFQ Intensity` = c(3100, 4100),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (identical(variant, "regular_intensity")) {
    names(data)[4:7] <- sub(" MaxLFQ Intensity$", " Intensity", names(data)[4:7])
  }
  if (identical(variant, "missing_intensity")) {
    data <- data[, c("Protein", "Protein ID", "Gene"), drop = FALSE]
  }
  if (identical(variant, "missing_protein")) {
    data <- data[, setdiff(names(data), "Protein ID"), drop = FALSE]
    data$Protein <- NULL
  }

  module_ci_prot_write_tsv(data, path)
}

module_ci_prot_run_values <- function(imported) {
  unique(as.character(imported$data[[imported$column_mapping$run_col]]))
}

module_ci_prot_feature_values <- function(imported) {
  unique(as.character(imported$data[[imported$column_mapping$protein_col]]))
}

module_ci_prot_assert_import_contract <- function(
    imported,
    expected_data_type,
    expected_samples,
    expected_features = NULL,
    expected_quantity_col
) {
  testthat::expect_identical(imported$data_type, expected_data_type)
  testthat::expect_identical(imported$column_mapping$quantity_col, expected_quantity_col)
  module_ci_assert_sample_identity(expected_samples, module_ci_prot_run_values(imported))
  if (!is.null(expected_features)) {
    module_ci_assert_feature_identity(expected_features, module_ci_prot_feature_values(imported))
  }
  invisible(imported)
}

module_ci_prot_workflow_data <- function() {
  state_manager <- new.env(parent = emptyenv())
  state_manager$workflow_type <- NULL
  state_manager$setWorkflowType <- function(type) {
    state_manager$workflow_type <- type
    invisible(type)
  }

  list2env(list(
    data_tbl = NULL,
    data_cln = NULL,
    data_format = NULL,
    data_type = NULL,
    column_mapping = NULL,
    fasta_file_path = NULL,
    aa_seq_tbl_final = NULL,
    fasta_metadata = NULL,
    taxon_id = NULL,
    organism_name = NULL,
    mixed_species_analysis = NULL,
    processing_log = list(),
    state_manager = state_manager,
    tab_status = list(
      setup_import = "pending",
      design_matrix = "disabled",
      quality_control = "disabled",
      normalization = "disabled",
      differential_expression = "disabled",
      enrichment_analysis = "disabled",
      session_summary = "disabled"
    )
  ), parent = emptyenv())
}

module_ci_prot_local_data <- function() {
  list2env(list(processing = TRUE, waiting_for_organism_selection = FALSE), parent = emptyenv())
}

module_ci_prot_expected_workflow_type <- function(format) {
  switch(format, pd_tmt = "TMT", diann = "DIA", "LFQ")
}

module_ci_prot_expected_report_template <- function(format) {
  switch(
    module_ci_prot_expected_workflow_type(format),
    DIA = "DIANN_report.rmd",
    TMT = "TMT_report.rmd",
    LFQ = "LFQ_report.rmd"
  )
}

module_ci_prot_assert_workflow_state <- function(workflow_data, format, expected_template = NULL) {
  expected_workflow_type <- module_ci_prot_expected_workflow_type(format)
  if (is.null(expected_template)) {
    expected_template <- module_ci_prot_expected_report_template(format)
  }

  testthat::expect_identical(workflow_data$state_manager$workflow_type, expected_workflow_type)
  testthat::expect_identical(workflow_data$data_format, format)
  testthat::expect_identical(workflow_data$tab_status$setup_import, "complete")
  testthat::expect_identical(workflow_data$tab_status$design_matrix, "pending")
  testthat::expect_identical(expected_template, module_ci_prot_expected_report_template(format))

  digest <- collect_state_digest(
    values = list(selected_omics = "proteomics"),
    workflow_states = list(proteomics = workflow_data)
  )
  testthat::expect_identical(digest$workflow_type_per_omic$proteomics, expected_workflow_type)
  testthat::expect_true("setup_import" %in% names(digest$export_paths$proteomics))
  invisible(digest)
}

module_ci_prot_apply_import_state <- function(imported, format, fasta_path = module_ci_prot_make_fasta()) {
  workflow_data <- module_ci_prot_workflow_data()
  local_data <- module_ci_prot_local_data()
  artifact_dir <- tempfile("module-ci-prot-import-artifacts-")
  dir.create(file.path(artifact_dir, "results"), recursive = TRUE)
  dir.create(file.path(artifact_dir, "scripts"), recursive = TRUE)

  checkpoints <- character()
  recordProtImportImportedData(
    dataImportResult = imported,
    captureCheckpoint = function(object, checkpoint_id, checkpoint_label) {
      checkpoints <<- c(checkpoints, paste(checkpoint_id, checkpoint_label, sep = ":"))
      invisible(object)
    },
    logInfo = function(...) NULL,
    messageFn = function(...) NULL
  )

  applyProtImportWorkflowWithStatus(
    workflowData = workflow_data,
    dataImportResult = imported,
    format = format,
    fastaPath = fasta_path,
    organismName = "Homo sapiens",
    experimentPaths = list(
      results_dir = file.path(artifact_dir, "results"),
      source_dir = file.path(artifact_dir, "scripts")
    ),
    sanitizeNames = FALSE,
    updateStatus = function(...) NULL,
    processFastaData = function(workflowData, fastaPath, organismName, experimentPaths, ...) {
      workflowData$aa_seq_tbl_final <- data.frame(
        uniprot_acc = c("P001", "P002"),
        sequence = c("MPEPTIDESEQ", "MPEPTIDESEQAA")
      )
      workflowData$fasta_metadata <- list(fasta_format = "standard", num_sequences = 2L)
      invisible(list(success = !is.null(fastaPath), cacheDir = file.path(experimentPaths$results_dir, "cache")))
    },
    logInfo = function(...) NULL,
    logWarn = function(...) NULL,
    showNotification = function(...) NULL
  )

  finalizeProtImportSetupState(
    workflowData = workflow_data,
    dataImportResult = imported,
    format = format,
    searchFilename = paste0(format, "_fixture.tsv"),
    fastaFilename = if (is.null(fasta_path)) NA_character_ else basename(fasta_path),
    taxonId = 9606L,
    organismName = "Homo sapiens",
    now = function() as.POSIXct("2026-05-05 00:00:00", tz = "UTC")
  )

  completeProtImportSuccessState(
    workflowData = workflow_data,
    localData = local_data,
    dataImportResult = imported,
    format = format,
    removeModal = function() NULL,
    showNotification = function(...) NULL,
    messageFn = function(...) NULL
  )

  list(workflow_data = workflow_data, checkpoints = checkpoints, artifact_dir = artifact_dir)
}
