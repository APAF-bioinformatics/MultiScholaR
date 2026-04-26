# fidelity-coverage-compare: shared
library(testthat)

makeFunctionWithOverrides <- function(fun, replacements) {
  fun_override <- fun
  environment(fun_override) <- list2env(replacements, parent = environment(fun))
  fun_override
}

restoreGlobalBinding <- function(name, had_value, old_value = NULL) {
  if (had_value) {
    assign(name, old_value, envir = .GlobalEnv)
  } else if (exists(name, envir = .GlobalEnv, inherits = FALSE)) {
    rm(list = name, envir = .GlobalEnv)
  }
}

test_that("createStudyParametersFile covers global fallbacks and missing contrast-column reporting", {
  config_env <- new.env(parent = emptyenv())
  config_env$config_list <- list(
    nested = list(alpha = 1, beta = c("x", "y")),
    internal_workflow_source_dir = "/tmp/skip-me"
  )

  source_dir <- tempfile("general-filemgmt-study-direct-")
  dir.create(source_dir, recursive = TRUE)
  withr::defer(unlink(source_dir, recursive = TRUE, force = TRUE))

  had_organism <- exists("organism_name", envir = .GlobalEnv, inherits = FALSE)
  old_organism <- if (had_organism) get("organism_name", envir = .GlobalEnv) else NULL
  had_taxon <- exists("taxon_id", envir = .GlobalEnv, inherits = FALSE)
  old_taxon <- if (had_taxon) get("taxon_id", envir = .GlobalEnv) else NULL
  withr::defer(restoreGlobalBinding("organism_name", had_organism, old_organism))
  withr::defer(restoreGlobalBinding("taxon_id", had_taxon, old_taxon))

  assign("organism_name", "Homo sapiens", envir = .GlobalEnv)
  assign("taxon_id", "9606", envir = .GlobalEnv)

  study_fun <- makeFunctionWithOverrides(
    createStudyParametersFile,
    list(
      requireNamespace = function(package, quietly = TRUE) FALSE
    )
  )

  study_file <- suppressMessages(study_fun(
    workflow_name = "proteomics",
    description = "direct study summary",
    source_dir_path = source_dir,
    contrasts_tbl = data.frame(name = "A-B", stringsAsFactors = FALSE),
    config_list_name = "config_list",
    env = config_env
  ))

  study_lines <- readLines(study_file)

  expect_identical(basename(study_file), "study_parameters.txt")
  expect_true(any(grepl("^Git Information:$", study_lines)))
  expect_true(any(grepl("Organism Name: Homo sapiens", study_lines, fixed = TRUE)))
  expect_true(any(grepl("Taxon ID: 9606", study_lines, fixed = TRUE)))
  expect_true(any(grepl("^Configuration Parameters:$", study_lines)))
  expect_true(any(grepl("Alpha: 1", study_lines, fixed = TRUE)))
  expect_true(any(grepl("\\[Column 'contrasts' not found in contrasts_tbl\\]", study_lines)))
  expect_false(any(grepl("internal_workflow_source_dir", study_lines, fixed = TRUE)))
})

test_that("createWorkflowArgsFromConfig records rich proteomics workflow metadata", {
  source_dir <- tempfile("general-filemgmt-workflow-direct-")
  dir.create(source_dir, recursive = TRUE)
  withr::defer(unlink(source_dir, recursive = TRUE, force = TRUE))

  had_config <- exists("config_list", envir = .GlobalEnv, inherits = FALSE)
  old_config <- if (had_config) get("config_list", envir = .GlobalEnv) else NULL
  had_organism <- exists("organism_name", envir = .GlobalEnv, inherits = FALSE)
  old_organism <- if (had_organism) get("organism_name", envir = .GlobalEnv) else NULL
  had_taxon <- exists("taxon_id", envir = .GlobalEnv, inherits = FALSE)
  old_taxon <- if (had_taxon) get("taxon_id", envir = .GlobalEnv) else NULL
  withr::defer(restoreGlobalBinding("config_list", had_config, old_config))
  withr::defer(restoreGlobalBinding("organism_name", had_organism, old_organism))
  withr::defer(restoreGlobalBinding("taxon_id", had_taxon, old_taxon))

  assign(
    "config_list",
    list(
      extra_section = list(mode = "strict"),
      internal_workflow_source_dir = "/tmp/skip-me"
    ),
    envir = .GlobalEnv
  )
  assign("organism_name", "Mus musculus", envir = .GlobalEnv)
  assign("taxon_id", "10090", envir = .GlobalEnv)

  final_object <- methods::new("ProteinQuantitativeData")
  final_object@args <- list(
    ruvIII_C_Varying = list(ruv_grouping_variable = "condition"),
    summary_section = list(
      flag = TRUE,
      many_flags = rep(c(TRUE, FALSE), 30),
      numeric_values = 1:6,
      labels = letters[1:6],
      frame = data.frame(x = 1:2),
      nothing = NULL
    )
  )

  workflow_file <- suppressMessages(suppressWarnings(createWorkflowArgsFromConfig(
    workflow_name = "proteomics",
    description = "direct workflow summary",
    source_dir_path = source_dir,
    final_s4_object = final_object,
    contrasts_tbl = data.frame(contrasts = "A-B", stringsAsFactors = FALSE),
    workflow_data = list(
      ruv_optimization_result = list(
        best_percentage = 15,
        best_k = 3,
        best_separation_score = 0.8123,
        best_composite_score = 0.7654,
        best_control_genes_index = c(TRUE, FALSE, TRUE, TRUE),
        separation_metric_used = "pearson",
        k_penalty_weight = 1.5,
        adaptive_k_penalty_used = TRUE,
        sample_size = 12
      ),
      fasta_metadata = list(
        fasta_format = "UniProt",
        num_sequences = 123,
        has_protein_evidence = TRUE,
        has_gene_names = TRUE,
        has_isoform_info = TRUE,
        has_status_info = FALSE
      ),
      mixed_species_analysis = list(
        enabled = TRUE,
        selected_organism = "Mouse",
        selected_taxon_id = "10090",
        filter_applied_at_import = TRUE,
        organism_distribution = data.frame(
          organism_name = c("Mouse", "Human"),
          taxon_id = c("10090", "9606"),
          protein_count = c(90, 10),
          percentage = c(90, 10),
          stringsAsFactors = FALSE
        )
      ),
      enrichment_organism_filter = list(
        enabled = TRUE,
        filter_applied = TRUE,
        target_taxon_id = "10090",
        proteins_before = 120,
        proteins_after = 96,
        proteins_removed = 24
      ),
      accession_cleanup_results = list(
        cleanup_applied = TRUE,
        aggregation_method = "best-score",
        delimiter_used = ":",
        proteins_before = 120,
        proteins_after = 96,
        had_full_metadata = TRUE
      ),
      protein_counts = list(
        after_qc_filtering = 88,
        after_ruv_filtering = 80,
        final_for_de = 72
      ),
      da_ui_params = list(
        q_value_threshold = 0.05,
        log_fold_change_cutoff = 1,
        treat_enabled = TRUE
      ),
      enrichment_ui_params = list(
        up_log2fc_cutoff = 1.2,
        down_log2fc_cutoff = -1.2,
        q_value_cutoff = 0.05,
        organism_selected = "Mouse",
        database_source = "GO"
      )
    )
  )))

  workflow_lines <- readLines(workflow_file)

  expect_identical(basename(workflow_file), "study_parameters.txt")
  expect_true(any(grepl("^Version Information:$", workflow_lines)))
  expect_true(any(grepl("^Organism Information:$", workflow_lines)))
  expect_true(any(grepl("Organism Name: Mus musculus", workflow_lines, fixed = TRUE)))
  expect_true(any(grepl("^Automatic RUV Optimization Results:$", workflow_lines)))
  expect_true(any(grepl("Best percentage: 15.0%", workflow_lines, fixed = TRUE)))
  expect_true(any(grepl("RUV grouping variable: condition", workflow_lines, fixed = TRUE)))
  expect_true(any(grepl("^FASTA File Processing:$", workflow_lines)))
  expect_true(any(grepl("^Mixed Species FASTA Analysis:$", workflow_lines)))
  expect_true(any(grepl("Mouse \\(Taxon 10090\\): 90 proteins", workflow_lines)))
  expect_true(any(grepl("^Enrichment Analysis - Organism Filtering:$", workflow_lines)))
  expect_true(any(grepl("Retention Rate: 80%", workflow_lines, fixed = TRUE)))
  expect_true(any(grepl("^Protein Accession Cleanup:$", workflow_lines)))
  expect_true(any(grepl("^Protein Filtering Summary:$", workflow_lines)))
  expect_true(any(grepl("^\\[summary_section\\]$", workflow_lines)))
  expect_true(any(grepl(
    "many_flags = logical vector \\[30 TRUE, 30 FALSE out of 60 total\\]",
    workflow_lines
  )))
  expect_true(any(grepl("numeric_values = numeric vector \\[6 values: 1, 2, 3, ...\\]", workflow_lines)))
  expect_true(any(grepl("frame = \\[Data frame: 2 rows x 1 cols - omitted for brevity\\]", workflow_lines)))
  expect_true(any(grepl("^  User Interface Parameters:$", workflow_lines)))
  expect_true(any(grepl("^Additional Configuration Parameters:$", workflow_lines)))
  expect_true(any(grepl("^Contrasts:$", workflow_lines)))
  expect_false(any(grepl("internal_workflow_source_dir", workflow_lines, fixed = TRUE)))
})

test_that("createWorkflowArgsFromConfig covers the RUV-skipped and no-config fallbacks", {
  source_dir <- tempfile("general-filemgmt-workflow-skip-")
  dir.create(source_dir, recursive = TRUE)
  withr::defer(unlink(source_dir, recursive = TRUE, force = TRUE))

  had_config <- exists("config_list", envir = .GlobalEnv, inherits = FALSE)
  old_config <- if (had_config) get("config_list", envir = .GlobalEnv) else NULL
  withr::defer(restoreGlobalBinding("config_list", had_config, old_config))
  assign("config_list", list(), envir = .GlobalEnv)

  workflow_file <- suppressMessages(suppressWarnings(createWorkflowArgsFromConfig(
    workflow_name = "proteomics",
    source_dir_path = source_dir,
    final_s4_object = methods::new("ProteinQuantitativeData"),
    workflow_data = list(
      ruv_optimization_result = list(ruv_skipped = TRUE)
    )
  )))

  workflow_lines <- readLines(workflow_file)

  expect_true(any(grepl("^RUV-III Batch Correction:$", workflow_lines)))
  expect_true(any(grepl("Status: Not Applied", workflow_lines, fixed = TRUE)))
  expect_true(any(grepl("No additional configuration parameters available", workflow_lines, fixed = TRUE)))
})
