library(testthat)

test_that("MCI-007.1 protein-level TMT/LFQ bypass starts protein QC from imported protein S4 state", {
  for (workflow_type in c("TMT", "LFQ")) {
    workflow <- module_ci_prot_protein_level_workflow(workflow_type = workflow_type)

    result <- runProteinS4CreationStep(
      workflow,
      logInfoFn = function(...) invisible(NULL)
    )

    expect_s4_class(result$proteinObj, "ProteinQuantitativeData")
    expect_identical(result$proteinIdCol, "Protein.Ids")
    expect_identical(result$proteinCount, 4L)
    expect_identical(workflow$state_manager$getHistory(), "protein_s4_initial")
    expect_identical(
      workflow$state_manager$saved()[[1L]]$config_object,
      list(s4_class = "ProteinQuantitativeData", protein_id_column = "Protein.Ids")
    )
    expect_identical(result$proteinObj@args$globalParameters$workflow_type, workflow_type)
    module_ci_prot_assert_protein_alignment(result$proteinObj)
  }

  invalid_workflow <- module_ci_prot_protein_level_workflow(protein_col = "MissingProtein")
  expect_error(
    runProteinS4CreationStep(invalid_workflow, logInfoFn = function(...) invisible(NULL)),
    "Protein ID column MissingProtein not found",
    fixed = TRUE
  )
  expect_identical(invalid_workflow$state_manager$getHistory(), character(0))
})

test_that("MCI-007.2 DIA rollup matrix preserves IQ/limpa metadata and peptide-count guards", {
  peptide_s4 <- module_ci_prot_peptide_rollup_object()
  workflow <- module_ci_prot_protein_workflow(peptide_s4, initial_name = "imputed")
  dirs <- list(
    peptide_qc_dir = file.path(tempdir(), "mci007-peptide-qc"),
    protein_qc_dir = file.path(tempdir(), "mci007-protein-qc")
  )
  dir.create(dirs$peptide_qc_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(dirs$protein_qc_dir, recursive = TRUE, showWarnings = FALSE)

  iq_result <- runProteinIqRollupApplyStep(
    workflow,
    experimentPaths = dirs,
    processLongFormatFn = module_ci_prot_iq_process_stub,
    readTsvFn = function(file, .name_repair) readr::read_tsv(file, show_col_types = FALSE),
    captureCheckpointFn = function(...) invisible(NULL),
    showNotificationFn = function(...) invisible(NULL),
    logInfoFn = function(...) invisible(NULL),
    logWarnFn = function(...) invisible(NULL),
    sleepFn = function(...) invisible(NULL)
  )

  expect_s4_class(iq_result$proteinObj, "ProteinQuantitativeData")
  expect_identical(workflow$state_manager$getHistory(), "protein_s4_created")
  expect_setequal(colnames(iq_result$proteinObj@protein_quant_table), c("Protein.Ids", "A1", "A2", "B1"))
  expect_setequal(iq_result$proteinObj@design_matrix$Run, c("A1", "A2", "B1"))
  expect_match(iq_result$resultText, "Algorithm: MaxLFQ", fixed = TRUE)
  module_ci_prot_assert_protein_alignment(iq_result$proteinObj)

  limpa_workflow <- module_ci_prot_protein_workflow(peptide_s4, initial_name = "imputed")
  limpa_result <- runProteinLimpaRollupApplyStep(
    limpa_workflow,
    experimentPaths = dirs,
    requireNamespaceFn = function(...) FALSE,
    isTestModeFn = function() TRUE,
    captureCheckpointFn = function(...) invisible(NULL),
    logInfoFn = function(...) invisible(NULL)
  )

  expect_s4_class(limpa_result$proteinObj, "ProteinQuantitativeData")
  expect_true(limpa_result$proteinObj@args$globalParameters$use_limpa)
  expect_identical(
    limpa_result$proteinObj@args$globalParameters$report_template,
    "DIANN_limpa_report.rmd"
  )
  expect_identical(
    limpa_workflow$state_manager$saved()[[1L]]$config_object$rollup_method,
    "limpa_dpc_quant"
  )
  expect_match(limpa_result$resultText, "DIANN_limpa_report.rmd", fixed = TRUE)

  production_workflow <- module_ci_prot_protein_workflow(peptide_s4, initial_name = "imputed")
  expect_error(
    runProteinLimpaRollupApplyStep(
      production_workflow,
      experimentPaths = dirs,
      requireNamespaceFn = function(...) FALSE,
      isTestModeFn = function() FALSE,
      captureCheckpointFn = function(...) invisible(NULL),
      logInfoFn = function(...) invisible(NULL)
    ),
    "limpa package is required",
    fixed = TRUE
  )
  expect_identical(production_workflow$state_manager$getHistory(), character(0))

  count_args <- limpa_result$proteinObj@args$limpa_dpc_quant_results$peptide_counts_per_protein
  limpa_result$proteinObj@args$limpa_dpc_quant_results$peptide_counts_per_protein <- count_args |>
    dplyr::mutate(
      peptide_count = ifelse(Protein.Ids == "P_SINGLE", 1L, peptide_count),
      peptidoform_count = ifelse(Protein.Ids == "P_SINGLE", 1L, peptidoform_count)
    )
  guarded <- filterMinNumPeptidesPerProtein(
    limpa_result$proteinObj,
    num_peptides_per_protein_thresh = 2,
    num_peptidoforms_per_protein_thresh = 2,
    verbose = FALSE
  )
  expect_false("P_SINGLE" %in% guarded@protein_quant_table$Protein.Ids)
})

test_that("MCI-007.3 accession cleanup covers delimiters, contaminants, isoforms, duplicates, and missing genes", {
  oracle <- module_ci_prot_accession_oracle()

  ranked <- suppressMessages(suppressWarnings(capture.output(
    accession_map <- chooseBestProteinAccessionHelper(
      input_tbl = oracle$input,
      acc_detail_tab = oracle$metadata,
      accessions_column = Protein.Ids,
      row_id_column = "uniprot_acc",
      group_id = row_id,
      delim = ";"
    )
  )))
  expect_type(ranked, "character")

  expect_false(any(grepl("CON__", accession_map$uniprot_acc, fixed = TRUE)))
  expect_match(accession_map$uniprot_acc[accession_map$row_id == 1], "P_BEST")
  expect_identical(accession_map$uniprot_acc[accession_map$row_id == 3], "P_ISO")
  expect_identical(accession_map$gene_names[accession_map$row_id == 4], "NA")

  object <- module_ci_prot_protein_object(
    module_ci_prot_protein_table(
      proteins = c("P1;P2", "P1;P2", "P3"),
      values = rbind(c(10, 10, 10, 10), c(20, 20, 20, 20), c(30, 30, 30, 30))
    )
  )
  workflow <- module_ci_prot_protein_workflow(object)
  cleaned_object <- object
  cleaned_object@protein_quant_table <- data.frame(
    Protein.Ids = c("P2", "P3"),
    A1 = c(15, 30),
    A2 = c(15, 30),
    B1 = c(15, 30),
    B2 = c(15, 30),
    check.names = FALSE
  )

  cleanup_result <- runProteinAccessionCleanupStep(
    workflow,
    delimiter = ";",
    aggregationMethod = "mean",
    chooseBestProteinAccessionFn = function(...) cleaned_object,
    nowFn = function() as.POSIXct("2026-05-05", tz = "UTC"),
    logInfoFn = function(...) invisible(NULL),
    logWarnFn = function(...) invisible(NULL),
    saveRdsFn = function(...) invisible(NULL),
    existsFn = function(...) TRUE,
    getFn = function(...) data.frame(database_id = c("P1", "P2", "P3"))
  )
  expect_true(cleanup_result$cleanupApplied)
  expect_identical(workflow$state_manager$getHistory(), "protein_accession_cleaned")
  expect_identical(workflow$qc_params$protein_qc$accession_cleanup$proteins_before, 2L)
  expect_identical(workflow$qc_params$protein_qc$accession_cleanup$proteins_after, 2L)

  skipped_workflow <- module_ci_prot_protein_workflow(object)
  skipped <- runProteinAccessionCleanupStep(
    skipped_workflow,
    delimiter = ";",
    aggregationMethod = "mean",
    nowFn = function() as.POSIXct("2026-05-05", tz = "UTC"),
    logInfoFn = function(...) invisible(NULL),
    logWarnFn = function(...) invisible(NULL),
    saveRdsFn = function(...) invisible(NULL),
    existsFn = function(...) FALSE
  )
  expect_false(skipped$cleanupApplied)
  expect_identical(skipped_workflow$state_manager$getHistory(), "protein_accession_cleaned")
})

test_that("MCI-007.4 protein filtering and dedup cover thresholds, non-finite values, duplicate choices, and tie-safe columns", {
  object <- module_ci_prot_protein_object(
    module_ci_prot_protein_table(
      proteins = c("P_KEEP", "P_DROP", "P_INF", "P_NAN"),
      values = rbind(
        c(20, 21, 22, 23),
        c(NA, 1, NA, 1),
        c(Inf, Inf, Inf, Inf),
        c(NaN, NaN, NaN, NaN)
      )
    ),
    args = list(
      removeRowsWithMissingValuesPercent = list(
        groupwise_percentage_cutoff = 49,
        max_groups_percentage_cutoff = 49,
        proteins_intensity_cutoff_percentile = 1
      )
    )
  )

  filtered <- suppressMessages(removeRowsWithMissingValuesPercent(object))
  expect_identical(filtered@protein_quant_table$Protein.Ids, "P_KEEP")
  module_ci_prot_assert_protein_alignment(filtered)

  dedup_object <- module_ci_prot_protein_object(
    module_ci_prot_protein_table(
      proteins = c("P_DUP", "P_DUP", "P_KEEP"),
      samples = c("Alpha", "Beta"),
      values = rbind(c(10, 20), c(30, 40), c(50, 60))
    ),
    design = module_ci_prot_protein_design(
      samples = c("Alpha", "Beta"),
      groups = c("A", "B"),
      replicates = c("R1", "R1")
    )
  )
  dedup_workflow <- module_ci_prot_protein_workflow(dedup_object)
  deduped <- runProteinDuplicateRemovalStep(
    dedup_workflow,
    aggregationMethod = "median",
    logInfoFn = function(...) invisible(NULL)
  )

  expect_setequal(colnames(deduped$deduplicatedS4@protein_quant_table), c("Protein.Ids", "Alpha", "Beta"))
  expect_equal(
    deduped$deduplicatedS4@protein_quant_table$Alpha[deduped$deduplicatedS4@protein_quant_table$Protein.Ids == "P_DUP"],
    20
  )
  expect_identical(dedup_workflow$state_manager$saved()[[1L]]$config_object$duplicates_found, "P_DUP")
})

test_that("MCI-007.5 protein replicate filter covers balanced, unbalanced, single-replicate, and disagreement cases", {
  object <- module_ci_prot_protein_object(
    module_ci_prot_protein_table(
      proteins = c("P_BALANCED", "P_SINGLE", "P_DISAGREE", "P_UNBALANCED"),
      values = rbind(
        c(10, 11, 12, 13),
        c(10, NA, NA, NA),
        c(30, NA, 31, NA),
        c(40, 41, 42, NA)
      )
    )
  )

  filtered <- suppressMessages(removeProteinsWithOnlyOneReplicate(
    object,
    core_utilisation = NA,
    grouping_variable = "group"
  ))
  expect_true("P_BALANCED" %in% filtered@protein_quant_table$Protein.Ids)
  expect_false("P_SINGLE" %in% filtered@protein_quant_table$Protein.Ids)
  expect_false("P_DISAGREE" %in% filtered@protein_quant_table$Protein.Ids)
  expect_false("P_UNBALANCED" %in% filtered@protein_quant_table$Protein.Ids)

  workflow <- module_ci_prot_protein_workflow(object)
  output_dir <- file.path(tempdir(), "mci007-protein-replicate")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  result <- runProteinReplicateFilterApplyStep(
    workflow,
    experimentPaths = list(protein_qc_dir = output_dir, source_dir = output_dir),
    groupingVariable = "group",
    parallelCores = 1,
    writeTsvFn = function(x, file) {
      readr::write_tsv(x, file)
      invisible(file)
    },
    saveRdsFn = function(...) invisible(NULL),
    logInfoFn = function(...) invisible(NULL),
    logWarnFn = function(...) invisible(NULL),
    isTestModeFn = function() TRUE
  )
  expect_identical(workflow$state_manager$getHistory(), "protein_replicate_filtered")
  expect_true(file.exists(result$outputFile))
  expect_identical(workflow$qc_params$protein_qc$replicate_filter$grouping_variable, "group")
  expect_identical(workflow$protein_counts$after_qc_filtering, result$proteinCount)
})

test_that("MCI-007.6 protein QC artifacts, counts, metadata, and state history remain coherent", {
  workflow <- module_ci_prot_protein_workflow()
  now <- as.POSIXct("2026-05-05 00:00:00", tz = "UTC")

  cleanup <- runProteinAccessionCleanupStep(
    workflow,
    delimiter = ";",
    aggregationMethod = "mean",
    nowFn = function() now,
    logInfoFn = function(...) invisible(NULL),
    logWarnFn = function(...) invisible(NULL),
    saveRdsFn = function(...) invisible(NULL),
    existsFn = function(...) FALSE
  )
  intensity <- runProteinIntensityFilterApplyStep(
    workflow,
    useStrictMode = FALSE,
    minRepsPerGroup = 1,
    minGroups = 1,
    intensityCutoffPercentile = 1,
    updateConfigParameterFn = function(theObject, function_name, parameter_name, new_value) {
      if (is.null(theObject@args[[function_name]])) {
        theObject@args[[function_name]] <- list()
      }
      theObject@args[[function_name]][[parameter_name]] <- new_value
      theObject
    },
    updateMissingValueParametersFn = function(theObject, min_reps_per_group, min_groups) {
      if (is.null(theObject@args$removeRowsWithMissingValuesPercent)) {
        theObject@args$removeRowsWithMissingValuesPercent <- list()
      }
      theObject@args$removeRowsWithMissingValuesPercent$groupwise_percentage_cutoff <- 50
      theObject@args$removeRowsWithMissingValuesPercent$max_groups_percentage_cutoff <- 50
      theObject
    },
    logInfoFn = function(...) invisible(NULL),
    nowFn = function() now
  )
  duplicate <- runProteinDuplicateRemovalStep(
    workflow,
    aggregationMethod = "mean",
    logInfoFn = function(...) invisible(NULL)
  )
  output_dir <- file.path(tempdir(), "mci007-state-artifacts")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  replicate <- runProteinReplicateFilterApplyStep(
    workflow,
    experimentPaths = list(protein_qc_dir = output_dir),
    groupingVariable = "group",
    parallelCores = 1,
    writeTsvFn = function(x, file) {
      readr::write_tsv(x, file)
      invisible(file)
    },
    saveRdsFn = function(...) invisible(NULL),
    logInfoFn = function(...) invisible(NULL),
    logWarnFn = function(...) invisible(NULL),
    isTestModeFn = function() TRUE,
    nowFn = function() now
  )

  expect_identical(
    workflow$state_manager$getHistory(),
    c(
      "protein_accession_cleaned",
      "protein_intensity_filtered",
      "duplicates_removed",
      "protein_replicate_filtered"
    )
  )
  expect_equal(nrow(intensity$filteredS4@protein_quant_table), 4L)
  expect_equal(duplicate$deduplicatedS4@protein_quant_table$Protein.Ids |> unique() |> length(), 3L)
  expect_equal(replicate$proteinCount, nrow(replicate$filteredS4@protein_quant_table))
  expect_false(cleanup$cleanupApplied)
  expect_identical(
    workflow$qc_params$protein_qc$intensity_filter$proteins_intensity_cutoff_percentile,
    1
  )
  expect_identical(workflow$qc_params$protein_qc$replicate_filter$timestamp, now)

  density_plot <- plotDensity(workflow$state_manager$getState(), grouping_variable = "group", title = "density")
  pca_plot <- plotPca(workflow$state_manager$getState(), grouping_variable = "group", label_column = NULL, title = "pca")
  rle_plot <- plotRle(workflow$state_manager$getState(), grouping_variable = "group")
  expect_s3_class(density_plot, "ggplot")
  expect_s3_class(pca_plot, "ggplot")
  expect_s3_class(rle_plot, "ggplot")
})
