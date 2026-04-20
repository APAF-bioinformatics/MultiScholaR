library(testthat)
library(shiny)

makeProtImportCharacterizationOverrides <- function(fun, replacements) {
  fun_override <- fun
  environment(fun_override) <- list2env(replacements, parent = environment(fun))
  fun_override
}

withProtImportCharacterizationGlobals <- function(object_names, code) {
  had_existing <- vapply(
    object_names,
    function(name) exists(name, envir = .GlobalEnv, inherits = FALSE),
    logical(1)
  )
  old_values <- lapply(seq_along(object_names), function(i) {
    if (had_existing[[i]]) {
      get(object_names[[i]], envir = .GlobalEnv, inherits = FALSE)
    } else {
      NULL
    }
  })
  names(old_values) <- object_names

  on.exit({
    for (name in rev(object_names)) {
      if (had_existing[[name]]) {
        assign(name, old_values[[name]], envir = .GlobalEnv)
      } else if (exists(name, envir = .GlobalEnv, inherits = FALSE)) {
        rm(list = name, envir = .GlobalEnv)
      }
    }
  }, add = TRUE)

  force(code)
}

makeProtImportCharacterizationFileInput <- function(path, type = "text/plain") {
  data.frame(
    name = basename(path),
    size = unname(file.info(path)$size),
    type = type,
    datapath = path,
    stringsAsFactors = FALSE
  )
}

makeProtImportCharacterizationFixture <- function(with_config = TRUE) {
  root_dir <- tempfile("prot-import-characterization-")
  dir.create(root_dir, recursive = TRUE)

  search_results_path <- file.path(root_dir, "report.tsv")
  writeLines(
    c(
      paste(
        c(
          "Protein.Group", "Protein.Ids", "Protein.Names", "Precursor.Id",
          "Modified.Sequence", "Stripped.Sequence", "Precursor.Charge",
          "Q.Value", "PG.Q.Value", "Run"
        ),
        collapse = "\t"
      ),
      paste(
        c(
          "P1", "P1", "Protein 1", "PEPTIDE1", "PEPTIDE1", "PEPTIDE1",
          "2", "0.01", "0.01", "Sample 1"
        ),
        collapse = "\t"
      )
    ),
    search_results_path
  )

  fasta_path <- file.path(root_dir, "proteins.fasta")
  writeLines(
    c(
      ">sp|P1|PROT1 OS=Homo sapiens OX=9606",
      "MPEPTIDESEQ"
    ),
    fasta_path
  )

  config_path <- NULL
  if (with_config) {
    config_path <- file.path(root_dir, "config.ini")
    writeLines("[generalParameters]\nmin_peptides_per_protein=2", config_path)
  }

  results_dir <- file.path(root_dir, "results")
  source_dir <- file.path(root_dir, "scripts")
  dir.create(results_dir, recursive = TRUE)
  dir.create(source_dir, recursive = TRUE)

  list(
    root_dir = root_dir,
    search_results_path = search_results_path,
    fasta_path = fasta_path,
    config_path = config_path,
    experiment_paths = list(
      results_dir = results_dir,
      source_dir = source_dir
    )
  )
}

makeProtImportCharacterizationWorkflow <- function() {
  shiny::reactiveValues(
    data_tbl = NULL,
    fasta_file_path = NULL,
    config_list = NULL,
    taxon_id = NULL,
    organism_name = NULL,
    design_matrix = NULL,
    data_cln = NULL,
    contrasts_tbl = NULL,
    state_manager = WorkflowState$new(),
    peptide_data = NULL,
    protein_log2_quant = NULL,
    protein_data = NULL,
    ruv_normalised_for_da_analysis_obj = NULL,
    da_analysis_results_list = NULL,
    uniprot_dat_cln = NULL,
    enrichment_results = NULL,
    data_format = NULL,
    data_type = NULL,
    column_mapping = NULL,
    aa_seq_tbl_final = NULL,
    fasta_metadata = NULL,
    uniprot_mapping = NULL,
    uniparc_mapping = NULL,
    mixed_species_analysis = NULL,
    tab_status = list(
      setup_import = "pending",
      design_matrix = "disabled",
      quality_control = "disabled",
      normalization = "disabled",
      differential_expression = "disabled",
      enrichment_analysis = "disabled",
      session_summary = "disabled"
    ),
    state_update_trigger = NULL,
    processing_log = list()
  )
}

mockProtImportCharacterizationResult <- function() {
  list(
    data = data.frame(
      Protein.Group = c("P1", "P2"),
      Protein.Ids = c("P1", "P2"),
      Precursor.Id = c("PEPTIDE1", "PEPTIDE2"),
      Run = c("Sample 1", "Sample 2"),
      Precursor.Quantity = c(1000, 2000),
      stringsAsFactors = FALSE
    ),
    data_type = "peptide",
    column_mapping = list(
      protein_col = "Protein.Group",
      peptide_col = "Precursor.Id",
      run_col = "Run",
      quantity_col = "Precursor.Quantity"
    )
  )
}

mockProtImportCharacterizationFasta <- function() {
  list(
    aa_seq_tbl_final = data.frame(
      uniprot_acc = c("P1", "P2"),
      sequence = c("MPEPTIDESEQ", "MPEPTIDESEQ2"),
      stringsAsFactors = FALSE
    ),
    fasta_metadata = list(
      fasta_format = "standard",
      num_sequences = 2
    )
  )
}

test_that("mod_prot_import_server preserves successful import side effects", {
  fixture <- makeProtImportCharacterizationFixture(with_config = TRUE)
  on.exit(unlink(fixture$root_dir, recursive = TRUE, force = TRUE), add = TRUE)

  workflow_data <- makeProtImportCharacterizationWorkflow()
  config_stub <- getDefaultProteomicsConfig()

  withProtImportCharacterizationGlobals(c("aa_seq_tbl_final", "config_list"), {
    server_under_test <- makeProtImportCharacterizationOverrides(
      mod_prot_import_server,
      list(
        requireNamespace = function(package, quietly = FALSE) {
          if (identical(package, "shinyFiles")) {
            FALSE
          } else {
            base::requireNamespace(package, quietly = quietly)
          }
        },
        importDIANNData = function(filepath, use_precursor_norm = TRUE) {
          expect_equal(filepath, fixture$search_results_path)
          expect_true(isTRUE(use_precursor_norm))
          mockProtImportCharacterizationResult()
        },
        processFastaFile = function(...) {
          mockProtImportCharacterizationFasta()
        },
        readConfigFile = function(file) {
          expect_equal(file, fixture$config_path)
          config_stub
        },
        .capture_checkpoint = function(...) invisible(NULL)
      )
    )

    testServer(
      server_under_test,
      args = list(
        workflow_data = workflow_data,
        experiment_paths = fixture$experiment_paths,
        volumes = NULL
      ),
      {
        session$setInputs(
          search_results_standard = makeProtImportCharacterizationFileInput(fixture$search_results_path),
          fasta_file_standard = makeProtImportCharacterizationFileInput(fixture$fasta_path),
          config_file_standard = makeProtImportCharacterizationFileInput(fixture$config_path),
          format_override = "diann",
          diann_use_precursor_norm = TRUE,
          sanitize_names = TRUE,
          mixed_species_fasta = FALSE,
          taxon_id = 9606,
          organism_name = "Homo sapiens",
          capture_checkpoints = FALSE
        )
        session$flushReact()

        session$setInputs(process_data = 1)
        session$flushReact()

        expect_equal(workflow_data$data_format, "diann")
        expect_equal(workflow_data$data_type, "peptide")
        expect_equal(workflow_data$state_manager$workflow_type, "DIA")
        expect_equal(workflow_data$tab_status$setup_import, "complete")
        expect_equal(workflow_data$fasta_file_path, fixture$fasta_path)
        expect_equal(sort(unique(workflow_data$data_tbl$Run)), c("sample_1", "sample_2"))
        expect_identical(workflow_data$data_cln, workflow_data$data_tbl)
        expect_equal(workflow_data$taxon_id, 9606)
        expect_equal(workflow_data$organism_name, "Homo sapiens")
        expect_equal(workflow_data$config_list, config_stub)
        expect_equal(workflow_data$fasta_metadata$num_sequences, 2)
        expect_false(workflow_data$mixed_species_analysis$enabled)
        expect_equal(workflow_data$processing_log$setup_import$detected_format, "diann")
        expect_equal(workflow_data$processing_log$setup_import$n_runs, 2)
        expect_true(exists("aa_seq_tbl_final", envir = .GlobalEnv, inherits = FALSE))
        expect_true(exists("config_list", envir = .GlobalEnv, inherits = FALSE))
      }
    )
  })
})

test_that("mod_prot_import_server preserves unsupported-format cleanup behavior", {
  fixture <- makeProtImportCharacterizationFixture(with_config = FALSE)
  on.exit(unlink(fixture$root_dir, recursive = TRUE, force = TRUE), add = TRUE)

  workflow_data <- makeProtImportCharacterizationWorkflow()
  workflow_data$data_tbl <- data.frame(old = 1)
  workflow_data$data_format <- "stale"
  workflow_data$data_type <- "stale"
  workflow_data$column_mapping <- list(old = TRUE)
  workflow_data$data_cln <- data.frame(old = 1)
  workflow_data$fasta_file_path <- "stale.fasta"
  workflow_data$aa_seq_tbl_final <- data.frame(old = 1)
  workflow_data$config_list <- list(old = TRUE)
  workflow_data$processing_log <- list(setup_import = list(old = TRUE))
  workflow_data$tab_status <- list(
    setup_import = "complete",
    design_matrix = "disabled",
    quality_control = "disabled",
    normalization = "disabled",
    differential_expression = "disabled",
    enrichment_analysis = "disabled",
    session_summary = "disabled"
  )

  server_under_test <- makeProtImportCharacterizationOverrides(
    mod_prot_import_server,
    list(
      requireNamespace = function(package, quietly = FALSE) {
        if (identical(package, "shinyFiles")) {
          FALSE
        } else {
          base::requireNamespace(package, quietly = quietly)
        }
      },
      .capture_checkpoint = function(...) invisible(NULL)
    )
  )

  testServer(
    server_under_test,
    args = list(
      workflow_data = workflow_data,
      experiment_paths = fixture$experiment_paths,
      volumes = NULL
    ),
    {
      session$setInputs(
        search_results_standard = makeProtImportCharacterizationFileInput(fixture$search_results_path),
        fasta_file_standard = makeProtImportCharacterizationFileInput(fixture$fasta_path),
        format_override = "unknown",
        mixed_species_fasta = FALSE,
        capture_checkpoints = FALSE
      )
      session$flushReact()

      session$setInputs(process_data = 1)
      session$flushReact()

      expect_null(workflow_data$data_tbl)
      expect_null(workflow_data$data_format)
      expect_null(workflow_data$data_type)
      expect_null(workflow_data$column_mapping)
      expect_null(workflow_data$data_cln)
      expect_null(workflow_data$fasta_file_path)
      expect_null(workflow_data$aa_seq_tbl_final)
      expect_null(workflow_data$config_list)
      expect_null(workflow_data$processing_log$setup_import)
      expect_equal(workflow_data$tab_status$setup_import, "incomplete")
    }
  )
})
