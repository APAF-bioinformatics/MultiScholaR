test_that("MCI-010.1 backend and taxon routing are explicit and validated", {
  supported <- buildProtEnrichSupportedOrganisms()

  gprofiler <- resolveProtEnrichAnalysisMethod(" 0009606 ", supported)
  expect_equal(gprofiler$analysisMethod, "gprofiler2")
  expect_true(isTRUE(gprofiler$organismSupported))
  expect_equal(gprofiler$methodInfo$species_id, "hsapiens")
  expect_equal(gprofiler$methodInfo$species_name, "Homo sapiens")

  mouse <- resolveProtEnrichAnalysisMethod("10090", supported)
  expect_equal(mouse$analysisMethod, "gprofiler2")
  expect_equal(mouse$methodInfo$species_id, "mmusculus")

  unsupported <- resolveProtEnrichAnalysisMethod("999999", supported)
  expect_equal(unsupported$analysisMethod, "clusterprofiler")
  expect_false(isTRUE(unsupported$organismSupported))
  expect_equal(unsupported$methodInfo$species_name, "Taxon ID 999999")

  expect_equal(normaliseProtEnrichTaxonId(9606), "9606")
  expect_error(normaliseProtEnrichTaxonId(""), "organism_taxid")
  expect_error(normaliseProtEnrichTaxonId("human"), "organism_taxid")
  expect_error(normaliseProtEnrichTaxonId("0"), "organism_taxid")

  input <- shiny::reactiveValues(organism_taxid = "9606")
  enrichment_data <- shiny::reactiveValues(
    analysis_method = NULL,
    organism_supported = NULL
  )
  current_method <- createProtEnrichCurrentAnalysisMethodReactive(
    input = input,
    enrichmentData = enrichment_data,
    supportedOrganismsFn = function() supported
  )

  expect_equal(shiny::isolate(current_method())$method, "gprofiler2")
  expect_equal(shiny::isolate(enrichment_data$analysis_method), "gprofiler2")
  expect_true(isTRUE(shiny::isolate(enrichment_data$organism_supported)))

  shiny::isolate(input$organism_taxid <- "999999")
  expect_equal(shiny::isolate(current_method())$method, "clusterprofiler")
  expect_equal(shiny::isolate(enrichment_data$analysis_method), "clusterprofiler")
  expect_false(isTRUE(shiny::isolate(enrichment_data$organism_supported)))
})

test_that("MCI-010.2 identifier and annotation matrices preserve controlled matching behavior", {
  uniprot_da <- module_ci_prot_enrich_da_results("uniprot")
  uniprot_object <- module_ci_prot_enrich_object(
    protein_id_column = "uniprot_acc",
    da_table = uniprot_da@da_data[[1]]
  )
  cluster_input <- resolveProtEnrichAnalysisInputColumns(
    methodInfo = module_ci_prot_enrich_method("clusterprofiler"),
    daResultsForEnrichment = uniprot_da,
    currentS4Object = uniprot_object
  )
  expect_equal(cluster_input$idColumn, "uniprot_acc")
  expect_equal(cluster_input$source, "s4_object")

  gprofiler_input <- resolveProtEnrichAnalysisInputColumns(
    methodInfo = module_ci_prot_enrich_method("gprofiler2"),
    daResultsForEnrichment = uniprot_da,
    currentS4Object = uniprot_object
  )
  expect_equal(gprofiler_input$idColumn, "gene_name")
  expect_true(isTRUE(gprofiler_input$geneNameOverrideApplied))

  ensembl_da <- module_ci_prot_enrich_da_results("ensembl_like")
  ensembl_da@da_data[[1]]$gene_name <- NULL
  ensembl_object <- module_ci_prot_enrich_object(
    protein_id_column = "ensembl_id",
    da_table = module_ci_prot_enrich_da_table("ensembl_like")
  )
  ensembl_input <- resolveProtEnrichAnalysisInputColumns(
    methodInfo = module_ci_prot_enrich_method("gprofiler2"),
    daResultsForEnrichment = ensembl_da,
    currentS4Object = ensembl_object
  )
  expect_equal(ensembl_input$idColumn, "ensembl_id")
  expect_false(isTRUE(ensembl_input$geneNameOverrideApplied))

  match_double <- function(da_results_s4,
                           uniprot_annotations,
                           protein_id_column,
                           uniprot_id_column,
                           gene_names_column) {
    ids <- toupper(as.character(da_results_s4@da_data[[1]][[protein_id_column]]))
    ids <- ids[!is.na(ids) & nzchar(ids)]
    annotation_ids <- toupper(as.character(uniprot_annotations[[uniprot_id_column]]))
    list(
      match_statistics = list(
        matched = sum(unique(ids) %in% annotation_ids),
        total = length(unique(ids)),
        match_rate = round(100 * sum(unique(ids) %in% annotation_ids) / max(length(unique(ids)), 1))
      ),
      gene_names_column = gene_names_column
    )
  }

  cases <- list(
    uniprot = 100,
    duplicate_ids = 100,
    missing_annotation = 100,
    mixed_case = 100,
    no_mappable_ids = 0
  )

  for (case_name in names(cases)) {
    da_results <- module_ci_prot_enrich_da_results(case_name)
    object <- module_ci_prot_enrich_object(da_table = da_results@da_data[[1]])
    annotations <- if (identical(case_name, "no_mappable_ids")) {
      module_ci_prot_enrich_annotations("uniprot")
    } else {
      module_ci_prot_enrich_annotations(case_name)
    }

    matched <- resolveProtEnrichAnnotationMatching(
      uniprotDatCln = annotations,
      daResultsForEnrichment = da_results,
      currentS4Object = object,
      matchAnnotationsFn = match_double,
      catFn = function(...) invisible(NULL)
    )

    expect_true(isTRUE(matched$attempted), info = case_name)
    expect_equal(matched$proteinIdColumn, "uniprot_acc", info = case_name)
    expect_equal(matched$matchRate, cases[[case_name]], info = case_name)
    expect_null(matched$warning, info = case_name)
  }
})

test_that("MCI-010.3 thresholds, filters, top-N, and empty-result boundaries are stable", {
  validated <- validateProtEnrichProcessParameters(
    organismTaxid = " 9606 ",
    upCutoff = "1",
    downCutoff = 0,
    qCutoff = "0.05",
    correctionMethod = " gSCS "
  )
  expect_identical(validated$organism_taxid, "9606")
  expect_equal(validated$up_cutoff, 1)
  expect_equal(validated$down_cutoff, 0)
  expect_equal(validated$q_cutoff, 0.05)
  expect_equal(validated$correction_method, "gSCS")

  expect_error(validateProtEnrichProcessParameters("9606", 1, 0, -0.01), "q_cutoff")
  expect_error(validateProtEnrichProcessParameters("9606", 1, 0, 1.01), "q_cutoff")
  expect_error(validateProtEnrichProcessParameters("9606", Inf, 0, 0.05), "up_cutoff")
  expect_error(validateProtEnrichProcessParameters("9606", 1, NA, 0.05), "down_cutoff")
  expect_error(validateProtEnrichProcessParameters("9606", 1, 0, 0.05, ""), "correction_method")

  args <- buildProtEnrichProcessEnrichmentsArgs(
    daResultsForEnrichment = module_ci_prot_enrich_da_results(),
    organismTaxid = "0009606",
    upCutoff = "1",
    downCutoff = "0",
    qCutoff = "0.05",
    pathwayDir = "/tmp/pathway",
    goAnnotations = module_ci_prot_enrich_annotations(),
    proteinIdColumn = "gene_name",
    contrastNames = "groupB-groupA",
    correctionMethod = "fdr"
  )
  expect_equal(args$processArgs$taxon_id, 9606)
  expect_equal(args$processArgs$up_cutoff, 1)
  expect_equal(args$processArgs$down_cutoff, 0)
  expect_equal(args$processArgs$q_cutoff, 0.05)
  expect_equal(args$processArgs$correction_method, "fdr")

  gprofiler_results <- rbind(
    module_ci_prot_enrich_gprofiler_table("positive"),
    module_ci_prot_enrich_gprofiler_table("negative")
  )
  gprofiler_table <- renderProtEnrichGprofilerResultsTable(
    gprofilerResults = gprofiler_results,
    directionFilter = "up",
    renderDtFn = function(expr) force(expr),
    datatableFn = module_ci_prot_enrich_datatable_capture,
    formatRoundFn = module_ci_prot_enrich_format_round_capture,
    catFn = function(...) invisible(NULL)
  )
  expect_true(all(gprofiler_table$data$directionality == "positive"))
  expect_true(all(c("pvalue", "p.adjust", "qvalue") %in% attr(gprofiler_table, "rounded_columns")))
  expect_equal(formatProtEnrichGprofilerSummaryText(gprofiler_results, "down"), "Showing 5 down-regulated pathways\n\nResults displayed in table below.")
  expect_equal(formatProtEnrichGprofilerSummaryText(gprofiler_results[0, ], "all"), "No gprofiler2 results available.")

  cluster_results <- rbind(
    module_ci_prot_enrich_cluster_table("up"),
    module_ci_prot_enrich_cluster_table("down")
  )
  cluster_table <- renderProtEnrichClusterProfilerResultsTable(
    clusterprofilerResults = cluster_results,
    directionFilter = "down",
    renderDtFn = function(expr) force(expr),
    datatableFn = module_ci_prot_enrich_datatable_capture,
    formatRoundFn = module_ci_prot_enrich_format_round_capture,
    catFn = function(...) invisible(NULL)
  )
  expect_true(all(cluster_table$data$directionality == "down"))
  expect_equal(formatProtEnrichClusterProfilerSummaryText(cluster_results, "up"), "Showing 3 up-regulated GO terms\n\nResults displayed in table below.")
  expect_equal(formatProtEnrichClusterProfilerSummaryText(cluster_results[0, ], "all"), "No clusterProfileR results available.")

  stringdb_table <- renderProtEnrichStringDbResultsTable(
    stringdbResults = data.frame(
      term = paste0("network_", 1:4),
      p_value = c(0.001, 0.02, 0.2, 0.03),
      stringsAsFactors = FALSE
    ),
    filterSignificant = TRUE,
    enrichmentPValThresh = 0.05,
    maxResults = 2,
    renderDtFn = function(expr) force(expr),
    datatableFn = module_ci_prot_enrich_datatable_capture
  )
  expect_equal(nrow(stringdb_table$data), 2)
  expect_true(all(stringdb_table$data$p_value < 0.05))

  source_filtered <- fake_gost(sources = "KEGG")$result
  expect_equal(nrow(source_filtered), 1)
  expect_equal(source_filtered$source, "KEGG")
})

test_that("MCI-010.4 deterministic backend doubles remain consumable", {
  da_results <- module_ci_prot_enrich_da_results()
  gprofiler_dir <- tempfile("module-ci-gprofiler-")
  cluster_dir <- tempfile("module-ci-cluster-")

  gprofiler_results <- runProtEnrichDeterministicProcessEnrichments(
    da_results = da_results,
    taxon_id = 9606,
    pathway_dir = gprofiler_dir,
    contrast_names = names(da_results@da_data)
  )
  expect_s4_class(gprofiler_results, "EnrichmentResults")
  expect_true(file.exists(file.path(gprofiler_dir, "groupB-groupA_up_enrichment_results.tsv")))
  expect_type(gprofiler_results@enrichment_data[["groupB-groupA"]]$up, "list")
  expect_true("result" %in% names(gprofiler_results@enrichment_data[["groupB-groupA"]]$up))

  gprofiler_collated <- buildProtEnrichAllContrastResults(
    enrichmentResults = gprofiler_results,
    methodInfo = module_ci_prot_enrich_method("gprofiler2"),
    catFn = function(...) invisible(NULL)
  )
  expect_true(all(gprofiler_collated[["groupB-groupA"]]$gprofiler_results$analysis_method == "gprofiler2"))
  expect_null(gprofiler_collated[["groupB-groupA"]]$clusterprofiler_results)

  cluster_results <- runProtEnrichDeterministicProcessEnrichments(
    da_results = da_results,
    taxon_id = 999999,
    pathway_dir = cluster_dir,
    contrast_names = names(da_results@da_data)
  )
  expect_s4_class(cluster_results, "EnrichmentResults")
  expect_true(file.exists(file.path(cluster_dir, "groupB-groupA_down_enrichment_results.tsv")))
  expect_s4_class(cluster_results@enrichment_data[["groupB-groupA"]]$up, "ProtEnrichDeterministicResult")

  cluster_collated <- buildProtEnrichAllContrastResults(
    enrichmentResults = cluster_results,
    methodInfo = module_ci_prot_enrich_method("clusterprofiler"),
    catFn = function(...) invisible(NULL)
  )
  expect_true(all(cluster_collated[["groupB-groupA"]]$clusterprofiler_results$analysis_method == "clusterprofiler"))
  expect_null(cluster_collated[["groupB-groupA"]]$gprofiler_results)

  expect_identical(
    resolveProtEnrichProcessEnrichmentsFn(
      processEnrichmentsFn = processEnrichments,
      deterministicProcessFn = runProtEnrichDeterministicProcessEnrichments,
      isTestModeFn = function() TRUE
    ),
    runProtEnrichDeterministicProcessEnrichments
  )
})

test_that("MCI-010.5 export ZIP and publication/report metadata preserve backend state", {
  gprofiler_archive <- module_ci_prot_enrich_archive("gprofiler2")
  expect_true("gprofiler2_results.tsv" %in% gprofiler_archive$listing$Name)
  expect_false("clusterProfileR_results.tsv" %in% gprofiler_archive$listing$Name)
  expect_true("enrichment_analysis_summary.txt" %in% gprofiler_archive$listing$Name)

  cluster_archive <- module_ci_prot_enrich_archive("clusterprofiler")
  expect_true("clusterProfileR_results.tsv" %in% cluster_archive$listing$Name)
  expect_false("gprofiler2_results.tsv" %in% cluster_archive$listing$Name)
  expect_true("enrichment_analysis_summary.txt" %in% cluster_archive$listing$Name)

  empty_archive <- module_ci_prot_enrich_archive("gprofiler2", result_case = "empty")
  expect_false("gprofiler2_results.tsv" %in% empty_archive$listing$Name)
  expect_true("enrichment_analysis_summary.txt" %in% empty_archive$listing$Name)

  extract_dir <- tempfile("module-ci-prot-enrich-archive-")
  dir.create(extract_dir)
  utils::unzip(cluster_archive$archive, exdir = extract_dir)
  summary_text <- paste(readLines(file.path(extract_dir, "enrichment_analysis_summary.txt"), warn = FALSE), collapse = "\n")
  expect_match(summary_text, "Contrast: B_vs_A", fixed = TRUE)
  expect_match(summary_text, "Analysis Method: clusterprofiler", fixed = TRUE)
  expect_match(summary_text, "Taxonomy ID: 999999", fixed = TRUE)
  expect_match(summary_text, "Q-value Cutoff: 0.05", fixed = TRUE)

  payload <- buildProtEnrichAnalysisResultsPayload(
    gprofilerResults = module_ci_prot_enrich_gprofiler_table("positive"),
    clusterprofilerResults = NULL,
    stringdbResults = NULL,
    fullEnrichmentResults = module_ci_prot_enrich_results("gprofiler2"),
    selectedContrast = "B_vs_A",
    analysisMethod = "gprofiler2",
    upCutoff = 1,
    downCutoff = 0,
    qCutoff = 0.05,
    organismTaxid = "9606"
  )
  expect_equal(payload$analysis_method, "gprofiler2")
  expect_equal(payload$selected_contrast, "B_vs_A")
  expect_equal(payload$parameters$organism_taxid, "9606")
  expect_equal(payload$parameters$q_cutoff, 0.05)

  skip_if_not_installed("openxlsx")
  paths <- module_ci_prot_enrich_paths()
  module_ci_prot_enrich_write_pathway_tables(paths, "gprofiler2", "B_vs_A")
  assign("module_ci_project_dirs", list(proteomics_MCI010 = paths, proteomics = paths), envir = .GlobalEnv)
  on.exit(rm("module_ci_project_dirs", envir = .GlobalEnv), add = TRUE)

  copyToResultsSummary(
    omic_type = "proteomics",
    experiment_label = "MCI010",
    contrasts_tbl = module_ci_prot_enrich_contrasts(),
    design_matrix = module_ci_prot_enrich_design(),
    force = TRUE,
    project_dirs_object_name = "module_ci_project_dirs"
  )

  workbook <- file.path(
    paths$results_summary_dir,
    "Publication_tables",
    "Pathway_enrichment_results_proteomics.xlsx"
  )
  expect_true(file.exists(workbook))
  index <- openxlsx::read.xlsx(workbook, sheet = "Enrichment_Index")
  index <- index[grepl("^Enrich_Sheet", index$Sheet), , drop = FALSE]
  expect_equal(index$Contrast, c("B_vs_A", "B_vs_A"))
  expect_setequal(index$Direction, c("up", "down"))
  sheet_data <- openxlsx::read.xlsx(workbook, sheet = index$Sheet[[1]])
  expect_equal(unique(sheet_data$analysis_method), "gprofiler2")
  expect_equal(unique(sheet_data$organism_taxid), 9606)
  expect_equal(unique(sheet_data$selected_contrast), "B_vs_A")
})

test_that("MCI-010.6 pre-seeded DA output initializes enrichment without rerunning DA", {
  workflow_data <- module_ci_prot_enrich_workflow()
  enrichment_data <- new.env(parent = emptyenv())
  enrichment_data$current_s4_object <- NULL
  enrichment_data$da_results_data <- NULL
  enrichment_data$contrasts_available <- NULL
  global_env <- new.env(parent = emptyenv())
  global_env$contrasts_tbl <- module_ci_prot_enrich_contrasts(
    raw = "groupB-groupA",
    friendly = "B_vs_A"
  )
  recorded <- new.env(parent = emptyenv())
  recorded$selected_update <- NULL
  recorded$notifications <- list()

  initialized <- registerProtEnrichSelectedTabObserver(
    selectedTabFn = function() "enrichment",
    workflowData = workflow_data,
    enrichmentData = enrichment_data,
    session = "module-ci-session",
    observeEventFn = function(eventExpr, handlerExpr, ignoreInit = FALSE) {
      eval(substitute(handlerExpr), parent.frame())
    },
    globalEnv = global_env,
    updateSelectInputFn = function(session, inputId, choices) {
      recorded$selected_update <- list(
        session = session,
        inputId = inputId,
        choices = choices
      )
    },
    showNotificationFn = function(...) {
      recorded$notifications <- c(recorded$notifications, list(list(...)))
    },
    catFn = function(...) invisible(NULL)
  )

  expect_true(isTRUE(initialized$initialized))
  expect_equal(initialized$reason, "initialized")
  expect_equal(initialized$source, "state_manager")
  expect_s4_class(enrichment_data$current_s4_object, "ProteinQuantitativeData")
  expect_identical(enrichment_data$da_results_data, workflow_data$da_analysis_results_list)
  expect_equal(enrichment_data$contrasts_available, "B_vs_A")
  expect_identical(recorded$selected_update$choices, c(B_vs_A = "B_vs_A"))
  expect_length(recorded$notifications, 0)

  input <- list(selected_contrast = "B_vs_A")
  preflight <- runProtEnrichObserverPreflight(
    input = input,
    enrichmentData = enrichment_data,
    workflowData = workflow_data,
    experimentPaths = module_ci_prot_enrich_paths(),
    currentAnalysisMethodFn = function() module_ci_prot_enrich_method("gprofiler2"),
    reqFn = function(selectedContrast, daResultsData) {
      expect_equal(selectedContrast, "B_vs_A")
      expect_identical(daResultsData, workflow_data$da_analysis_results_list)
      invisible(TRUE)
    },
    handoffObserverRunFn = function(input,
                                    enrichmentData,
                                    workflowData,
                                    experimentPaths,
                                    currentAnalysisMethodFn) {
      list(
        selectedContrast = input$selected_contrast,
        method = currentAnalysisMethodFn()$method,
        source = "preseeded_da"
      )
    }
  )

  expect_equal(preflight$selectedContrast, "B_vs_A")
  expect_equal(preflight$handoffResult$method, "gprofiler2")
  expect_equal(preflight$handoffResult$source, "preseeded_da")
})
