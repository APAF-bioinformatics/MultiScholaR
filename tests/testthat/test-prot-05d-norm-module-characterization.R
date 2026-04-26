# fidelity-coverage-compare: shared
library(testthat)
library(shiny)

if (!methods::isClass("mockProtNormCharacterizationData")) {
  methods::setClass(
    "mockProtNormCharacterizationData",
    slots = c(
      protein_quant_table = "data.frame",
      protein_id_column = "character",
      args = "list"
    )
  )
}

makeProtNormCharacterizationData <- function(label = "prot_norm_fixture") {
  methods::new(
    "mockProtNormCharacterizationData",
    protein_quant_table = data.frame(
      Protein.Ids = c("p1", "p2"),
      S1 = c(10, 20),
      S2 = c(30, 40),
      stringsAsFactors = FALSE
    ),
    protein_id_column = "Protein.Ids",
    args = list(
      label = label,
      globalParameters = list(workflow_type = "DIA")
    )
  )
}

multiScholaRNamespace <- function() {
  asNamespace("MultiScholaR")
}

hasMultiScholaRBinding <- function(name) {
  exists(name, envir = multiScholaRNamespace(), inherits = FALSE)
}

getMultiScholaRBinding <- function(name) {
  get(name, envir = multiScholaRNamespace(), inherits = FALSE)
}

localMultiScholaRBinding <- function(name, value, .local_envir = parent.frame()) {
  target_env <- multiScholaRNamespace()
  had_binding <- exists(name, envir = target_env, inherits = FALSE)
  old_value <- if (had_binding) get(name, envir = target_env, inherits = FALSE) else NULL
  was_locked <- had_binding && bindingIsLocked(name, target_env)

  if (was_locked) {
    unlockBinding(name, target_env)
  }
  assign(name, value, envir = target_env)
  if (was_locked) {
    lockBinding(name, target_env)
  }

  withr::defer({
    if (exists(name, envir = target_env, inherits = FALSE) &&
        bindingIsLocked(name, target_env)) {
      unlockBinding(name, target_env)
    }
    if (had_binding) {
      assign(name, old_value, envir = target_env)
    } else if (exists(name, envir = target_env, inherits = FALSE)) {
      rm(list = name, envir = target_env)
    }
    if (was_locked && exists(name, envir = target_env, inherits = FALSE)) {
      lockBinding(name, target_env)
    }
  }, envir = .local_envir)
}

localMultiScholaRBindings <- function(bindings, .local_envir = parent.frame()) {
  for (name in names(bindings)) {
    localMultiScholaRBinding(name, bindings[[name]], .local_envir = .local_envir)
  }

  invisible(NULL)
}

skipIfMissingMultiScholaRBindings <- function(...) {
  names <- unlist(list(...), use.names = FALSE)
  missing <- names[!vapply(names, hasMultiScholaRBinding, logical(1))]
  if (length(missing) > 0) {
    skip(paste("requires extracted helper bindings:", paste(missing, collapse = ", ")))
  }
}

makeProtNormCharacterizationHarness <- function(
  current_state = "ruv_corrected",
  current_s4 = makeProtNormCharacterizationData()
) {
  capture <- new.env(parent = emptyenv())
  capture$saved_states <- list()
  capture$reverted_states <- character()

  root_dir <- tempfile("prot-norm-characterization-")
  dir.create(root_dir, recursive = TRUE)
  protein_qc_dir <- file.path(root_dir, "qc")
  source_dir <- file.path(root_dir, "source")
  dir.create(protein_qc_dir, recursive = TRUE)
  dir.create(source_dir, recursive = TRUE)

  state_manager <- new.env(parent = emptyenv())
  state_manager$current_state <- current_state
  state_manager$getState <- function(state) current_s4
  state_manager$saveState <- function(...) {
    capture$saved_states[[length(capture$saved_states) + 1]] <<- list(...)
  }
  state_manager$getHistory <- function() {
    c("raw_data_s4", "sample_filtered", "protein_replicate_filtered", current_state)
  }
  state_manager$revertToState <- function(state) {
    capture$reverted_states <<- c(capture$reverted_states, state)
    current_s4
  }

  workflow_data <- new.env(parent = emptyenv())
  workflow_data$design_matrix <- NULL
  workflow_data$state_manager <- state_manager
  workflow_data$tab_status <- list(normalization = "pending", differential_expression = "disabled")
  workflow_data$state_update_trigger <- NULL
  workflow_data$protein_counts <- list()
  workflow_data$ruv_optimization_result <- list(ruv_skipped = FALSE)
  workflow_data$contrasts_tbl <- data.frame(friendly_names = "A vs B", stringsAsFactors = FALSE)
  workflow_data$config_list <- list(globalParameters = list(workflow_type = "DIA"))
  workflow_data$fasta_metadata <- list(database = "ref")
  workflow_data$accession_cleanup_results <- list(clean = TRUE)
  workflow_data$qc_params <- list(minimum = 1)
  workflow_data$mixed_species_analysis <- NULL

  experiment_paths <- list(
    protein_qc_dir = protein_qc_dir,
    source_dir = source_dir
  )

  list(
    workflow_data = workflow_data,
    experiment_paths = experiment_paths,
    session = shiny::MockShinySession$new(),
    current_s4 = current_s4,
    capture = capture
  )
}

launchProtNormModule <- function(harness, selected_tab = NULL) {
  shiny::withReactiveDomain(harness$session, {
    mod_prot_norm_server(
      id = "norm",
      workflow_data = harness$workflow_data,
      experiment_paths = harness$experiment_paths,
      omic_type = "proteomics",
      experiment_label = "Proteomics",
      selected_tab = selected_tab
    )
  })

  harness$session$getReturned()
}

setProtNormInputs <- function(session, ...) {
  values <- list(...)
  names(values) <- paste0("norm-", names(values))
  do.call(session$setInputs, values)
  shiny:::flushReact()
}

test_that("ProteinQuantitativeData chooseBestProteinAccession covers accession cleanup", {
  protein_object <- methods::new(
    "ProteinQuantitativeData",
    protein_quant_table = data.frame(
      Protein.Ids = c("P1;ALT", "P1"),
      S1 = c(1, 2),
      S2 = c(3, 4),
      stringsAsFactors = FALSE
    ),
    design_matrix = data.frame(
      Run = c("S1", "S2"),
      group = c("A", "B"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    group_id = "group",
    protein_id_column = "Protein.Ids",
    protein_id_table = data.frame(
      Protein.Ids = c("P1;ALT", "P1"),
      gene_names = c("GENE1", "GENE1"),
      stringsAsFactors = FALSE
    ),
    args = list(chooseBestProteinAccession = list())
  )

  method_env <- environment(methods::selectMethod("chooseBestProteinAccession", "ProteinQuantitativeData")@.Data)
  local_mocked_bindings(
    chooseBestProteinAccessionHelper = function(input_tbl,
                                                acc_detail_tab,
                                                accessions_column,
                                                row_id_column,
                                                group_id,
                                                delim) {
      data.frame(
        row_id = input_tbl$row_id,
        uniprot_acc = "P1",
        stringsAsFactors = FALSE
      )
    },
    rankProteinAccessionHelper = function(input_tbl,
                                          acc_detail_tab,
                                          accessions_column,
                                          row_id_column,
                                          group_id,
                                          delim) {
      data.frame(
        Protein.Ids = "P1",
        seqinr_accession_column = "P1;ALT",
        num_gene_names = 1,
        gene_names = "GENE1",
        is_unique = TRUE,
        stringsAsFactors = FALSE
      )
    },
    .env = method_env
  )

  accession_method <- methods::selectMethod("chooseBestProteinAccession", "ProteinQuantitativeData")
  cleaned <- accession_method(
    protein_object,
    delim = ";",
    seqinr_obj = list(aa_seq_tbl_final = data.frame(uniprot_acc = "P1", stringsAsFactors = FALSE)),
    seqinr_accession_column = "uniprot_acc",
    replace_zero_with_na = FALSE,
    aggregation_method = "sum"
  )

  expect_s4_class(cleaned, "ProteinQuantitativeData")
  expect_equal(cleaned@protein_quant_table$Protein.Ids, "P1")
  expect_equal(cleaned@protein_quant_table$S1, 3)
  expect_equal(cleaned@protein_quant_table$S2, 7)
  expect_true("Protein.Ids_list" %in% names(cleaned@protein_id_table))
})

test_that("createGridQC GridPlotData method assembles plot sections and save path", {
  grid <- methods::new(
    "GridPlotData",
    pca_plots = list(),
    density_plots = list(),
    rle_plots = list(),
    pearson_plots = list(),
    cancor_plots = list(),
    limpa_plots = list(),
    pca_titles = list(),
    density_titles = list(),
    rle_titles = list(),
    pearson_titles = list(),
    cancor_titles = list(),
    limpa_titles = list()
  )
  base_plot <- ggplot2::ggplot(data.frame(x = 1:2, y = 2:1), ggplot2::aes(x, y)) +
    ggplot2::geom_point()

  grid@pca_plots <- list(
    pca_plot_before_cyclic_loess_group = base_plot,
    pca_plot_before_ruvIIIc_group = base_plot,
    pca_plot_after_ruvIIIc_group = base_plot
  )
  grid@density_plots <- list(density_plot_before_cyclic_loess_group = base_plot)
  grid@rle_plots <- list(rle_plot_before_cyclic_loess_group = base_plot)
  grid@pearson_plots <- list(pearson_correlation_pair_before_cyclic_loess = base_plot)
  grid@cancor_plots <- list(cancor_plot_before_ruvIIIc = base_plot)
  grid@limpa_plots <- list(
    dpc_curve = base_plot,
    missing_comparison = base_plot,
    intensity_distribution = base_plot,
    summary = base_plot
  )
  grid@pca_titles <- list("PCA before", "PCA RUV", "PCA after")
  grid@density_titles <- list("Density before")
  grid@rle_titles <- list("RLE before")
  grid@pearson_titles <- list("Pearson before")
  grid@cancor_titles <- list("Cancor before")
  grid@limpa_titles <- list("Limpa summary")

  captured <- new.env(parent = emptyenv())
  method_env <- environment(methods::selectMethod("createGridQC", "GridPlotData")@.Data)
  local_mocked_bindings(
    ggsave = function(plot, filename, ...) {
      dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
      writeLines("mock plot", filename)
      captured$filename <- filename
      invisible(filename)
    },
    .env = method_env
  )

  output_dir <- tempfile("grid-qc-")
  result <- methods::selectMethod("createGridQC", "GridPlotData")(
    grid,
    ncol = 2,
    save_path = output_dir,
    file_name = "qc-grid"
  )

  expect_true(inherits(result, "patchwork") || inherits(result, "ggplot"))
  expect_true(file.exists(file.path(output_dir, "qc-grid.png")))
  expect_identical(captured$filename, file.path(output_dir, "qc-grid.png"))
})

test_that("mod_prot_norm_server preserves apply-correlation public behavior", {
  harness <- makeProtNormCharacterizationHarness()
  final_s4 <- makeProtNormCharacterizationData(label = "final")

  localMultiScholaRBindings(
    list(
      pearsonCorForSamplePairs = function(...) {
        data.frame(sample1 = "S1", sample2 = "S2", pearson = 0.99, stringsAsFactors = FALSE)
      },
      filterSamplesByProteinCorrelationThreshold = function(...) final_s4,
      updateProteinFiltering = function(...) "filter_plot",
      plotPca = function(...) "pca_plot"
    )
  )

  state <- launchProtNormModule(harness)
  state$ruv_normalized_obj <- harness$current_s4
  state$ruv_optimization_result <- list(ruv_skipped = FALSE)

  setProtNormInputs(
    harness$session,
    ruv_mode = "automatic",
    ruv_grouping_variable = "group",
    min_pearson_correlation_threshold = 0.8
  )
  setProtNormInputs(harness$session, apply_correlation_filter = 1)

  expect_true(shiny::isolate(state$correlation_filtering_complete))
  expect_identical(shiny::isolate(state$correlation_filtered_obj), final_s4)
  expect_identical(harness$workflow_data$ruv_normalised_for_da_analysis_obj, final_s4)
  expect_equal(harness$workflow_data$tab_status$normalization, "complete")
  expect_equal(harness$workflow_data$tab_status$differential_expression, "pending")
  expect_length(harness$capture$saved_states, 1)
  expect_equal(harness$capture$saved_states[[1]]$state_name, "correlation_filtered")
})

test_that("mod_prot_norm_server preserves skip-correlation public behavior", {
  harness <- makeProtNormCharacterizationHarness()

  localMultiScholaRBindings(
    list(
      updateProteinFiltering = function(...) "filter_plot",
      plotPca = function(...) "pca_plot"
    )
  )

  state <- launchProtNormModule(harness)
  state$ruv_normalized_obj <- harness$current_s4
  state$ruv_optimization_result <- list(ruv_skipped = FALSE)

  setProtNormInputs(
    harness$session,
    ruv_mode = "automatic",
    ruv_grouping_variable = "group"
  )
  setProtNormInputs(harness$session, skip_correlation_filter = 1)

  expect_true(shiny::isolate(state$correlation_filtering_complete))
  expect_identical(shiny::isolate(state$correlation_filtered_obj), harness$current_s4)
  expect_identical(harness$workflow_data$ruv_normalised_for_da_analysis_obj, harness$current_s4)
  expect_equal(harness$workflow_data$tab_status$differential_expression, "pending")
  expect_length(harness$capture$saved_states, 1)
  expect_true(isTRUE(harness$capture$saved_states[[1]]$config_object$skipped))
})

test_that("mod_prot_norm_server preserves export prereq and success behavior", {
  harness <- makeProtNormCharacterizationHarness(current_state = "correlation_filtered")

  state <- launchProtNormModule(harness)

  setProtNormInputs(
    harness$session,
    norm_method = "cyclicloess",
    ruv_mode = "automatic",
    export_filtered_session = 1
  )

  expect_false(file.exists(file.path(harness$experiment_paths$source_dir, "filtered_session_data_latest.rds")))

  state$correlation_filtering_complete <- TRUE
  state$correlation_filtered_obj <- harness$current_s4
  state$best_k <- 2
  state$correlation_threshold <- 0.75

  setProtNormInputs(harness$session, export_filtered_session = 2)

  expect_true(file.exists(file.path(harness$experiment_paths$source_dir, "filtered_session_data_latest.rds")))
  expect_true(file.exists(file.path(harness$experiment_paths$source_dir, "filtered_session_summary.txt")))
})

test_that("mod_prot_norm_server preserves reset behavior", {
  harness <- makeProtNormCharacterizationHarness()

  state <- launchProtNormModule(harness)
  state$normalization_complete <- TRUE
  state$ruv_complete <- TRUE
  state$correlation_filtering_complete <- TRUE
  state$normalized_protein_obj <- harness$current_s4
  state$ruv_normalized_obj <- harness$current_s4
  state$correlation_filtered_obj <- harness$current_s4
  state$best_k <- 2
  state$control_genes_index <- c(TRUE, FALSE)
  state$correlation_vector <- data.frame(value = 0.95)
  state$correlation_threshold <- 0.8
  state$final_qc_plot <- "pca_plot"
  state$final_filtering_plot <- "filter_plot"
  state$post_norm_filtering_plot <- "post_plot"
  state$filtering_summary_text <- "summary"
  state$ruv_optimization_result <- list(best_k = 2)

  setProtNormInputs(
    harness$session,
    ruv_mode = "automatic",
    reset_normalization = 1
  )

  expect_false(shiny::isolate(state$normalization_complete))
  expect_false(shiny::isolate(state$ruv_complete))
  expect_false(shiny::isolate(state$correlation_filtering_complete))
  expect_null(shiny::isolate(state$normalized_protein_obj))
  expect_null(shiny::isolate(state$ruv_normalized_obj))
  expect_null(shiny::isolate(state$correlation_filtered_obj))
  expect_equal(harness$workflow_data$tab_status$normalization, "pending")
  expect_equal(harness$workflow_data$tab_status$differential_expression, "disabled")
  expect_equal(harness$capture$reverted_states[[1]], "protein_replicate_filtered")
})

test_that("mod_prot_norm_server preserves invalid-tab warning behavior", {
  harness <- makeProtNormCharacterizationHarness(current_state = "intensity_filtered")
  selected_tab_value <- shiny::reactiveVal("overview")

  state <- launchProtNormModule(
    harness,
    selected_tab = function() selected_tab_value()
  )

  selected_tab_value("normalization")
  shiny:::flushReact()

  expect_false(shiny::isolate(state$pre_norm_qc_generated))
  expect_equal(harness$workflow_data$tab_status$differential_expression, "disabled")
})

test_that("mod_prot_norm_server preserves pre-normalization error handling", {
  harness <- makeProtNormCharacterizationHarness(current_state = "protein_replicate_filtered")
  selected_tab_value <- shiny::reactiveVal("overview")

  localMultiScholaRBindings(
    list(
      plotPca = function(...) stop("pre-qc boom")
    )
  )

  state <- launchProtNormModule(
    harness,
    selected_tab = function() selected_tab_value()
  )

  selected_tab_value("normalization")
  shiny:::flushReact()

  expect_false(shiny::isolate(state$pre_norm_qc_generated))
  expect_false(file.exists(file.path(harness$experiment_paths$protein_qc_dir, "pre_norm_pca.png")))
})

test_that("mod_prot_norm_server preserves normalization error handling", {
  harness <- makeProtNormCharacterizationHarness(current_state = "protein_replicate_filtered")

  localMultiScholaRBindings(
    list(
      normaliseBetweenSamples = function(...) stop("normalization boom")
    )
  )

  state <- launchProtNormModule(harness)

  setProtNormInputs(
    harness$session,
    norm_method = "cyclicloess",
    ruv_mode = "automatic",
    run_normalization = 1
  )

  expect_false(shiny::isolate(state$normalization_complete))
  expect_false(shiny::isolate(state$ruv_complete))
  expect_equal(harness$workflow_data$tab_status$differential_expression, "disabled")
})

test_that("mod_prot_norm_server preserves correlation error handling", {
  harness <- makeProtNormCharacterizationHarness()

  state <- launchProtNormModule(harness)
  state$ruv_normalized_obj <- NULL

  setProtNormInputs(
    harness$session,
    ruv_mode = "automatic",
    ruv_grouping_variable = "group",
    min_pearson_correlation_threshold = 0.8,
    apply_correlation_filter = 1
  )

  expect_false(shiny::isolate(state$correlation_filtering_complete))
  expect_length(harness$capture$saved_states, 0)
})

test_that("mod_prot_norm_server preserves export error handling", {
  harness <- makeProtNormCharacterizationHarness(current_state = "correlation_filtered")
  harness$experiment_paths$source_dir <- file.path(harness$experiment_paths$source_dir, "missing")

  state <- launchProtNormModule(harness)
  state$correlation_filtering_complete <- TRUE
  state$correlation_filtered_obj <- harness$current_s4
  state$best_k <- 2
  state$correlation_threshold <- 0.75

  setProtNormInputs(
    harness$session,
    norm_method = "cyclicloess",
    ruv_mode = "automatic",
    export_filtered_session = 1
  )

  expect_false(file.exists(file.path(harness$experiment_paths$source_dir, "filtered_session_data_latest.rds")))
})

test_that("mod_prot_norm_server preserves reset error handling", {
  harness <- makeProtNormCharacterizationHarness()
  harness$workflow_data$state_manager$getHistory <- function() stop("reset boom")

  state <- launchProtNormModule(harness)
  state$normalization_complete <- TRUE
  state$ruv_complete <- TRUE

  setProtNormInputs(
    harness$session,
    ruv_mode = "automatic",
    reset_normalization = 1
  )

  expect_true(shiny::isolate(state$normalization_complete))
  expect_true(shiny::isolate(state$ruv_complete))
  expect_equal(harness$workflow_data$tab_status$differential_expression, "disabled")
})

test_that("prot norm shared helper builders preserve QC contracts", {
  skipIfMissingMultiScholaRBindings(
    "generateProtNormCompositeFromFiles",
    "buildProtNormDensityPlot",
    "buildProtNormCorrelationPlot",
    "resolveProtNormQcStateObject",
    "updateProtNormFinalFilteringPlot",
    "updateProtNormFinalQcPlot"
  )
  skip_if_not_installed("patchwork")
  skip_if_not_installed("png")
  skip_if_not_installed("ggplot2")

  generate_composite <- getMultiScholaRBinding("generateProtNormCompositeFromFiles")
  build_density <- getMultiScholaRBinding("buildProtNormDensityPlot")
  build_correlation <- getMultiScholaRBinding("buildProtNormCorrelationPlot")
  resolve_state <- getMultiScholaRBinding("resolveProtNormQcStateObject")
  update_filter_plot <- getMultiScholaRBinding("updateProtNormFinalFilteringPlot")
  update_qc_plot <- getMultiScholaRBinding("updateProtNormFinalQcPlot")

  img_path <- tempfile(fileext = ".png")
  png::writePNG(array(1, dim = c(1, 1, 4)), img_path)
  composite <- generate_composite(
    plotFiles = c(img_path),
    ncol = 1,
    rowLabels = list(pca = "a)"),
    columnLabels = "Pre"
  )
  expect_type(composite, "list")
  expect_true(all(c("plot", "width", "height") %in% names(composite)))

  pca_calls <- list()
  box_calls <- list()
  density_plot <- build_density(
    s4Object = "s4",
    aesthetics = list(color_var = "group", shape_var = "batch"),
    plotPcaFn = function(...) {
      pca_calls <<- list(...)
      "pca_plot"
    },
    plotPcaBoxFn = function(...) {
      box_calls <<- list(...)
      "density_plot"
    }
  )
  expect_equal(density_plot, "density_plot")
  expect_equal(pca_calls$grouping_variable, "group")
  expect_equal(box_calls$grouping_variable, "group")

  pearson_calls <- list()
  correlation_plot <- build_correlation(
    s4Object = "s4",
    colorVar = "batch",
    isRunningFn = function() TRUE,
    withProgressFn = function(message, detail, value, expr) force(expr),
    incProgressFn = function(value) invisible(value),
    plotPearsonFn = function(...) {
      pearson_calls <<- list(...)
      "corr_plot"
    }
  )
  expect_equal(correlation_plot, "corr_plot")
  expect_equal(pearson_calls$correlation_group, "batch")

  state_manager <- new.env(parent = emptyenv())
  state_manager$current_state <- "normalized"
  state_manager$getState <- function(state) {
    expect_equal(state, "normalized")
    "s4_obj"
  }
  resolved <- resolve_state(
    stateManager = state_manager,
    reqFn = function(x) x,
    messageFn = function(...) NULL
  )
  expect_equal(resolved$currentS4, "s4_obj")

  final_s4 <- makeProtNormCharacterizationData(label = "final")
  norm_data <- new.env(parent = emptyenv())
  filtering_plot <- update_filter_plot(
    finalS4ForDe = final_s4,
    normData = norm_data,
    omicType = "proteomics",
    experimentLabel = "exp1",
    updateProteinFilteringFn = function(data, step_name, omic_type, experiment_label, return_grid, overwrite) {
      expect_equal(step_name, "12_correlation_filtered")
      expect_equal(omic_type, "proteomics")
      expect_equal(experiment_label, "exp1")
      expect_true(return_grid)
      expect_true(overwrite)
      "filter_plot"
    },
    messageFn = function(...) NULL
  )
  qc_plot <- update_qc_plot(
    finalS4ForDe = final_s4,
    normData = norm_data,
    title = "Final Data",
    getPlotAestheticsFn = function() list(color_var = "group", shape_var = "batch"),
    plotPcaFn = function(object, grouping_variable, label_column, shape_variable, title, font_size) {
      expect_identical(object, final_s4)
      expect_equal(grouping_variable, "group")
      expect_equal(shape_variable, "batch")
      expect_equal(title, "Final Data")
      expect_equal(font_size, 8)
      "qc_plot"
    },
    messageFn = function(...) NULL
  )
  expect_equal(filtering_plot, "filter_plot")
  expect_equal(qc_plot, "qc_plot")
})

test_that("prot norm shared normalization workflow preserves orchestration", {
  skipIfMissingMultiScholaRBindings("runProtNormNormalizationWorkflow")
  run_normalization_workflow <- getMultiScholaRBinding("runProtNormNormalizationWorkflow")

  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- "state_manager"
  norm_data <- new.env(parent = emptyenv())
  calls <- list()
  progress <- character()
  notification <- NULL

  result <- run_normalization_workflow(
    input = list(norm_method = "cyclicloess", ruv_mode = "automatic"),
    workflowData = workflow_data,
    normData = norm_data,
    experimentPaths = list(protein_qc_dir = tempdir(), source_dir = tempdir()),
    omicType = "proteomics",
    experimentLabel = "exp1",
    checkMemoryUsageFn = function(threshold_gb, context) {
      calls$memory <<- list(threshold_gb = threshold_gb, context = context)
    },
    generatePostNormalizationQcFn = function(object) {
      calls$post_qc_generated <<- object
    },
    generateRuvCorrectedQcFn = function(object) {
      calls$ruv_qc_generated <<- object
    },
    getRuvGroupingVariableFn = function() "group",
    prepareNormalizationRunFn = function(stateManager, normData) {
      calls$prepare <<- list(stateManager = stateManager, normData = normData)
      list(currentS4 = "current_obj")
    },
    runBetweenSamplesStepFn = function(currentS4, normMethod, normData, proteinQcDir) {
      calls$between <<- list(currentS4 = currentS4, normMethod = normMethod, normData = normData, proteinQcDir = proteinQcDir)
      "normalized_obj"
    },
    runPostNormalizationQcStepFn = function(normalizedS4, normData, generatePostNormalizationQcFn) {
      calls$post_step <<- list(normalizedS4 = normalizedS4, normData = normData)
      generatePostNormalizationQcFn(normalizedS4)
    },
    resolveRuvParametersFn = function(normalizedS4, input, normData, workflowData, sourceDir, getRuvGroupingVariableFn) {
      calls$resolve_ruv <<- list(normalizedS4 = normalizedS4, sourceDir = sourceDir, grouping = getRuvGroupingVariableFn())
    },
    applyRuvCorrectionStepFn = function(normalizedS4, normData, getRuvGroupingVariableFn) {
      calls$apply_ruv <<- list(normalizedS4 = normalizedS4, grouping = getRuvGroupingVariableFn())
      "ruv_obj"
    },
    finalizeRuvCleanupStepFn = function(ruvCorrectedS4, input, normData, workflowData, omicType, experimentLabel) {
      calls$cleanup <<- list(ruvCorrectedS4 = ruvCorrectedS4, omicType = omicType, experimentLabel = experimentLabel)
      "ruv_clean"
    },
    resolveStep6QcObjectFn = function(step5Object, normData) {
      calls$resolve_step6 <<- list(step5Object = step5Object, normData = normData)
      "step6_obj"
    },
    runStep6RuvQcFn = function(ruvMode, step6Object, normData, qcDir, generateRuvCorrectedQcFn) {
      calls$step6_qc <<- list(ruvMode = ruvMode, step6Object = step6Object, qcDir = qcDir)
      generateRuvCorrectedQcFn(step6Object)
    },
    generateCompositeQcFigureFn = function(ruvMode, qcDir, omicType) {
      calls$composite <<- list(ruvMode = ruvMode, qcDir = qcDir, omicType = omicType)
    },
    finalizeWorkflowStateFn = function(normData) {
      calls$finalize <<- normData
    },
    buildCompletionNotificationFn = function(ruvMode) {
      calls$completion_message <<- ruvMode
      "done"
    },
    withProgressFn = function(message, value, expr) force(expr),
    incProgressFn = function(value, detail) {
      progress <<- c(progress, detail)
    },
    showNotificationFn = function(message, type, duration) {
      notification <<- list(message = message, type = type, duration = duration)
    },
    messageFn = function(...) NULL
  )

  expect_true(isTRUE(result))
  expect_equal(calls$memory$context, "Normalization Start")
  expect_identical(calls$prepare$stateManager, "state_manager")
  expect_equal(calls$between$normMethod, "cyclicloess")
  expect_equal(calls$post_qc_generated, "normalized_obj")
  expect_equal(calls$resolve_ruv$grouping, "group")
  expect_equal(calls$apply_ruv$grouping, "group")
  expect_equal(calls$cleanup$ruvCorrectedS4, "ruv_obj")
  expect_equal(calls$resolve_step6$step5Object, "ruv_clean")
  expect_equal(calls$step6_qc$ruvMode, "automatic")
  expect_equal(calls$ruv_qc_generated, "step6_obj")
  expect_equal(calls$composite$omicType, "proteomics")
  expect_identical(calls$finalize, norm_data)
  expect_equal(notification$message, "done")
  expect_equal(
    progress,
    c(
      "Normalizing between samples...",
      "Generating post-normalization QC plots...",
      "Determining RUV parameters...",
      "Applying RUV-III batch correction...",
      "Generating RUV-corrected QC plots..."
    )
  )
})

test_that("prot norm shared apply-correlation workflow preserves orchestration", {
  skipIfMissingMultiScholaRBindings("runProtNormApplyCorrelationWorkflow")
  run_apply_correlation_workflow <- getMultiScholaRBinding("runProtNormApplyCorrelationWorkflow")

  ruv_s4 <- makeProtNormCharacterizationData(label = "ruv")
  final_s4 <- makeProtNormCharacterizationData(label = "final")
  norm_data <- new.env(parent = emptyenv())
  norm_data$ruv_optimization_result <- list(ruv_skipped = FALSE)
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$ruv_optimization_result <- list(ruv_skipped = FALSE)
  progress_calls <- character()
  calls <- list()

  result <- run_apply_correlation_workflow(
    ruvS4 = ruv_s4,
    correlationThreshold = 0.7,
    normData = norm_data,
    experimentPaths = list(protein_qc_dir = tempdir()),
    workflowData = workflow_data,
    omicType = "proteomics",
    experimentLabel = "exp1",
    getRuvGroupingVariableFn = function() "group",
    getPlotAestheticsFn = function() list(color_var = "group", shape_var = "batch"),
    runVectorStepFn = function(...) {
      calls$vector <<- list(...)
      data.frame(sample = "S1", correlation = 0.9)
    },
    runFilterStepFn = function(...) {
      calls$filter <<- list(...)
      final_s4
    },
    updateFinalFilteringPlotFn = function(...) {
      calls$updateFiltering <<- list(...)
      "filter_plot"
    },
    updateFinalQcPlotFn = function(...) {
      calls$updateQc <<- list(...)
      "qc_plot"
    },
    saveCorrelationResultsFn = function(...) {
      calls$save <<- list(...)
      invisible(NULL)
    },
    withProgressFn = function(message, value, expr) force(expr),
    incProgressFn = function(value, detail) {
      progress_calls <<- c(progress_calls, detail)
    },
    gcFn = function() {
      calls$gc <<- TRUE
    },
    messageFn = function(...) NULL
  )

  expect_identical(result, final_s4)
  expect_equal(progress_calls, c(
    "Calculating sample correlations...",
    "Filtering low-correlation samples...",
    "Updating tracking...",
    "Saving results..."
  ))
  expect_equal(calls$vector$correlationThreshold, 0.7)
  expect_identical(calls$filter$correlationVec$correlation, 0.9)
  expect_identical(calls$updateFiltering$finalS4ForDe, final_s4)
  expect_identical(calls$updateQc$finalS4ForDe, final_s4)
  expect_identical(calls$save$finalS4ForDe, final_s4)
  expect_true(isTRUE(calls$gc))
})

test_that("prot norm shared export helpers preserve summary and artifact workflows", {
  skipIfMissingMultiScholaRBindings(
    "buildProtNormExportSummaryContent",
    "saveProtNormExportArtifacts",
    "runProtNormExportSessionWorkflow"
  )
  build_export_summary <- getMultiScholaRBinding("buildProtNormExportSummaryContent")
  save_export_artifacts <- getMultiScholaRBinding("saveProtNormExportArtifacts")
  run_export_workflow <- getMultiScholaRBinding("runProtNormExportSessionWorkflow")

  session_data <- list(
    final_protein_count = 12,
    final_sample_count = 4,
    contrasts_tbl = data.frame(friendly_names = c("A vs B", "A vs C"), stringsAsFactors = FALSE),
    normalization_method = "cyclicloess",
    ruv_mode = "automatic",
    ruv_k = 2,
    correlation_threshold = 0.75,
    accession_cleanup_results = list(clean = TRUE),
    qc_params = list(min = 1),
    fasta_metadata = list(db = "ref"),
    ruv_optimization_result = list(best_k = 2),
    mixed_species_analysis = list(enabled = TRUE)
  )

  summary_text <- build_export_summary(
    sessionData = session_data,
    sessionFilename = "filtered_session_data_20260411_120000.rds",
    timeFn = function() as.POSIXct("2026-04-11 12:00:00", tz = "UTC"),
    formatTimeFn = function(x, format) if (identical(format, "%Y-%m-%d %H:%M:%S")) "2026-04-11 12:00:00" else "20260411_120000"
  )
  expect_match(summary_text, "Proteins: 12", fixed = TRUE)

  saved_rds <- list()
  written_summary <- NULL
  artifacts <- save_export_artifacts(
    sessionData = session_data,
    sourceDir = tempdir(),
    timeFn = function() as.POSIXct("2026-04-11 12:00:00", tz = "UTC"),
    formatTimeFn = function(x, format) if (identical(format, "%Y%m%d_%H%M%S")) "20260411_120000" else "2026-04-11 12:00:00",
    saveRdsFn = function(object, path) {
      saved_rds[[length(saved_rds) + 1]] <<- path
    },
    writeLinesFn = function(text, con) {
      written_summary <<- list(text = text, con = con)
    },
    fileExistsFn = function(path) FALSE,
    messageFn = function(...) NULL
  )
  expect_equal(artifacts$summaryFilepath, file.path(tempdir(), "filtered_session_summary.txt"))
  expect_equal(written_summary$con, artifacts$summaryFilepath)
  expect_true(any(grepl("filtered_session_data_20260411_120000.rds", saved_rds, fixed = TRUE)))

  progress_calls <- character()
  workflow_result <- run_export_workflow(
    workflowData = list(
      contrasts_tbl = session_data$contrasts_tbl,
      config_list = list(globalParameters = list(workflow_type = "DIA")),
      accession_cleanup_results = session_data$accession_cleanup_results,
      qc_params = session_data$qc_params,
      fasta_metadata = session_data$fasta_metadata,
      ruv_optimization_result = session_data$ruv_optimization_result,
      mixed_species_analysis = session_data$mixed_species_analysis,
      state_manager = list(current_state = "correlation_filtered", getState = function(state) makeProtNormCharacterizationData(label = state))
    ),
    normData = list(best_k = 2, correlation_threshold = 0.75),
    sourceDir = tempdir(),
    input = list(norm_method = "cyclicloess", ruv_mode = "automatic"),
    collectSessionDataFn = function(...) session_data,
    saveExportArtifactsFn = function(...) {
      progress_calls <<- c(progress_calls, "save")
      artifacts
    },
    withProgressFn = function(message, value, expr) force(expr),
    incProgressFn = function(value, detail) progress_calls <<- c(progress_calls, detail),
    messageFn = function(...) NULL
  )
  expect_identical(workflow_result$sessionData, session_data)
  expect_identical(workflow_result$exportArtifacts, artifacts)
  expect_true("save" %in% progress_calls)
})

test_that("prot norm shared QC generation helpers preserve artifact recording", {
  skipIfMissingMultiScholaRBindings(
    "generateProtNormPreNormalizationQcArtifacts",
    "generateProtNormPostNormalizationQcArtifacts",
    "generateProtNormRuvCorrectedQcArtifacts",
    "regenerateProtNormQcForAestheticChange",
    "initializeProtNormQcPlotPaths",
    "recordProtNormQcPlotPath",
    "saveProtNormQcPlotArtifact"
  )
  generate_pre_qc <- getMultiScholaRBinding("generateProtNormPreNormalizationQcArtifacts")
  generate_post_qc <- getMultiScholaRBinding("generateProtNormPostNormalizationQcArtifacts")
  generate_ruv_qc <- getMultiScholaRBinding("generateProtNormRuvCorrectedQcArtifacts")
  regenerate_qc <- getMultiScholaRBinding("regenerateProtNormQcForAestheticChange")
  initialize_paths <- getMultiScholaRBinding("initializeProtNormQcPlotPaths")
  record_path <- getMultiScholaRBinding("recordProtNormQcPlotPath")
  save_artifact <- getMultiScholaRBinding("saveProtNormQcPlotArtifact")

  saved_files <- character()
  save_stub <- function(qcDir, filename, plotObject, width, height, dpi = 150, savePlotFn = NULL) {
    saved_files <<- c(saved_files, filename)
    file.path(qcDir, filename)
  }
  record_stub <- function(qcPlotPaths, stage, plotType, path) {
    qcPlotPaths <- initialize_paths(qcPlotPaths)
    qcPlotPaths[[stage]][[plotType]] <- path
    qcPlotPaths
  }

  expect_true("post_filtering" %in% names(initialize_paths(NULL)))
  paths <- record_path(NULL, "post_filtering", "pca", "file.png")
  expect_equal(paths$post_filtering$pca, "file.png")
  expect_null(save_artifact(NULL, "plot.png", "plot", 8, 6, savePlotFn = function(...) stop("should not save")))

  pre_result <- generate_pre_qc(
    stateManager = list(
      current_state = "protein_replicate_filtered",
      getState = function(state) makeProtNormCharacterizationData(label = state)
    ),
    qcDir = tempdir(),
    aesthetics = list(color_var = "group", shape_var = "batch"),
    qcPlotPaths = NULL,
    messageFn = function(...) NULL,
    gcFn = function(...) NULL,
    plotPcaFn = function(...) "pca",
    plotRleFn = function(...) "rle",
    buildDensityFn = function(...) "density",
    buildCorrelationFn = function(...) "correlation",
    saveArtifactFn = save_stub,
    recordPathFn = record_stub
  )
  expect_true(all(c("pca", "rle", "density", "correlation") %in% names(pre_result$post_filtering)))

  post_result <- generate_post_qc(
    normalizedS4 = makeProtNormCharacterizationData(label = "post"),
    qcDir = tempdir(),
    aesthetics = list(color_var = "group", shape_var = "batch"),
    qcPlotPaths = pre_result,
    messageFn = function(...) NULL,
    gcFn = function(...) NULL,
    plotPcaFn = function(...) "pca",
    plotRleFn = function(...) "rle",
    buildDensityFn = function(...) "density",
    buildCorrelationFn = function(...) "correlation",
    saveArtifactFn = save_stub,
    recordPathFn = record_stub
  )
  expect_true(all(c("pca", "rle", "density", "correlation") %in% names(post_result$post_normalization)))

  ruv_result <- generate_ruv_qc(
    ruvCorrectedS4 = makeProtNormCharacterizationData(label = "ruv"),
    qcDir = tempdir(),
    aesthetics = list(color_var = "group", shape_var = "batch"),
    qcPlotPaths = post_result,
    messageFn = function(...) NULL,
    gcFn = function(...) NULL,
    plotPcaFn = function(...) "pca",
    plotRleFn = function(...) "rle",
    buildDensityFn = function(...) "density",
    buildCorrelationFn = function(...) "correlation",
    saveArtifactFn = save_stub,
    recordPathFn = record_stub
  )
  expect_true(all(c("pca", "rle", "density", "correlation") %in% names(ruv_result$ruv_corrected)))

  calls <- character()
  regenerate_qc(
    normData = list(
      pre_norm_qc_generated = TRUE,
      normalization_complete = TRUE,
      normalized_protein_obj = "norm_obj",
      ruv_complete = TRUE,
      ruv_normalized_obj = "ruv_obj"
    ),
    generatePreNormalizationQcFn = function() calls <<- c(calls, "pre"),
    generatePostNormalizationQcFn = function(x) calls <<- c(calls, paste("post", x)),
    generateRuvCorrectedQcFn = function(x) calls <<- c(calls, paste("ruv", x)),
    messageFn = function(...) NULL
  )
  expect_true(all(c("pre", "post norm_obj", "ruv ruv_obj") %in% calls))
})

test_that("prot norm shared correlation persistence helpers preserve state handoff", {
  skipIfMissingMultiScholaRBindings(
    "resolveProtNormCorrelationResultFilenames",
    "saveProtNormCorrelationResults",
    "finalizeProtNormCorrelationWorkflowState",
    "runProtNormCorrelationVectorStep",
    "runProtNormCorrelationFilterStep",
    "prepareProtNormSkippedCorrelationState"
  )
  resolve_filenames <- getMultiScholaRBinding("resolveProtNormCorrelationResultFilenames")
  save_results <- getMultiScholaRBinding("saveProtNormCorrelationResults")
  finalize_state <- getMultiScholaRBinding("finalizeProtNormCorrelationWorkflowState")
  run_vector_step <- getMultiScholaRBinding("runProtNormCorrelationVectorStep")
  run_filter_step <- getMultiScholaRBinding("runProtNormCorrelationFilterStep")
  prepare_skipped_state <- getMultiScholaRBinding("prepareProtNormSkippedCorrelationState")

  skipped_files <- resolve_filenames(
    normRuvOptimizationResult = list(ruv_skipped = TRUE),
    workflowRuvOptimizationResult = NULL,
    messageFn = function(...) NULL
  )
  expect_true(skipped_files$ruvWasSkipped)

  final_s4 <- makeProtNormCharacterizationData(label = "correlation-final")
  saved_tsv <- NULL
  saved_rds <- NULL
  file_names <- save_results(
    finalS4ForDe = final_s4,
    experimentPaths = list(protein_qc_dir = tempdir()),
    normRuvOptimizationResult = list(ruv_skipped = TRUE),
    workflowRuvOptimizationResult = NULL,
    writeTsvFn = function(data, path) {
      saved_tsv <<- list(data = data, path = path)
    },
    saveRdsFn = function(object, path) {
      saved_rds <<- list(object = object, path = path)
    },
    messageFn = function(...) NULL
  )
  expect_equal(saved_tsv$data, final_s4@protein_quant_table)
  expect_identical(saved_rds$object, final_s4)
  expect_true(file_names$ruvWasSkipped)

  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- new.env(parent = emptyenv())
  saved_state <- NULL
  workflow_data$state_manager$saveState <- function(...) saved_state <<- list(...)
  workflow_data$tab_status <- list(normalization = "pending", differential_expression = "disabled")
  workflow_data$state_update_trigger <- as.POSIXct("2026-04-11 10:00:00", tz = "UTC")
  metrics <- finalize_state(
    finalS4ForDe = final_s4,
    workflowData = workflow_data,
    correlationThreshold = 0.65,
    skipped = FALSE,
    timeFn = function() as.POSIXct("2026-04-11 10:05:00", tz = "UTC"),
    messageFn = function(...) NULL,
    catFn = function(...) NULL
  )
  expect_equal(saved_state$state_name, "correlation_filtered")
  expect_equal(metrics$finalProteinCount, 2)

  norm_data <- new.env(parent = emptyenv())
  correlation_vec <- run_vector_step(
    ruvS4 = final_s4,
    correlationThreshold = 0.7,
    normData = norm_data,
    getRuvGroupingVariableFn = function() "group",
    pearsonCorForSamplePairsFn = function(object, tech_rep_remove_regex, correlation_group) {
      data.frame(sample = "S1", correlation = 0.8)
    },
    timeFn = function() as.POSIXct("2026-04-11 10:00:00", tz = "UTC"),
    messageFn = function(...) NULL
  )
  expect_equal(correlation_vec$correlation, 0.8)
  expect_equal(norm_data$correlation_threshold, 0.7)

  filtered_s4 <- run_filter_step(
    ruvS4 = final_s4,
    correlationVec = correlation_vec,
    correlationThreshold = 0.75,
    normData = norm_data,
    filterSamplesFn = function(ruvS4, pearson_correlation_per_pair, min_pearson_correlation_threshold) final_s4,
    gcFn = function() NULL,
    timeFn = function() as.POSIXct("2026-04-11 10:00:02", tz = "UTC"),
    messageFn = function(...) NULL
  )
  expect_identical(filtered_s4, final_s4)

  skipped_s4 <- prepare_skipped_state(
    ruvS4 = final_s4,
    normData = norm_data,
    gcFn = function() NULL,
    messageFn = function(...) NULL
  )
  expect_identical(skipped_s4, final_s4)
})

test_that("prot norm shared plotting and choice helpers preserve defaults and artifacts", {
  skipIfMissingMultiScholaRBindings(
    "prepareProtNormOptimizationResultsTable",
    "getProtNormRuvCanonicalCorrelationPlot",
    "updateProtNormDesignDrivenChoices",
    "getProtNormPlotAesthetics",
    "getProtNormRuvGroupingVariable",
    "buildProtNormLabelPlot",
    "buildProtNormTitlePlot",
    "loadProtNormImageAsPlot",
    "shouldProtNormAutoGeneratePreQc"
  )
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("png")

  prepare_table <- getMultiScholaRBinding("prepareProtNormOptimizationResultsTable")
  get_cancor_plot <- getMultiScholaRBinding("getProtNormRuvCanonicalCorrelationPlot")
  update_choices <- getMultiScholaRBinding("updateProtNormDesignDrivenChoices")
  get_plot_aesthetics <- getMultiScholaRBinding("getProtNormPlotAesthetics")
  get_ruv_grouping <- getMultiScholaRBinding("getProtNormRuvGroupingVariable")
  build_label_plot <- getMultiScholaRBinding("buildProtNormLabelPlot")
  build_title_plot <- getMultiScholaRBinding("buildProtNormTitlePlot")
  load_image_plot <- getMultiScholaRBinding("loadProtNormImageAsPlot")
  should_auto_pre_qc <- getMultiScholaRBinding("shouldProtNormAutoGeneratePreQc")

  optimization_table <- prepare_table(list(
    optimization_results = data.frame(
      percentage = c(5, 10),
      separation_score = c(0.12345, 0.56789),
      composite_score = c(0.24681, 0.97531)
    ),
    best_percentage = 10
  ))
  expect_true(isTRUE(optimization_table$hasResults))
  expect_equal(optimization_table$bestPercentage, 10)
  expect_equal(optimization_table$results$separation_score, c(0.1235, 0.5679))

  expect_equal(get_cancor_plot(list(best_cancor_plot = "plot_obj")), "plot_obj")
  expect_null(get_cancor_plot(list(best_cancor_plot = NULL)))
  expect_equal(get_plot_aesthetics(NULL, ""), list(color_var = "group", shape_var = "group"))
  expect_equal(get_ruv_grouping(NULL), "group")
  expect_equal(get_ruv_grouping("batch"), "batch")
  expect_true(should_auto_pre_qc("normalization", "protein_replicate_filtered", FALSE))
  expect_false(should_auto_pre_qc("qc", "protein_replicate_filtered", FALSE))

  design_calls <- list()
  design_messages <- character()
  update_choices(
    session = "session",
    designMatrix = data.frame(
      group = c("A", "B"),
      batch = c("x", "y"),
      sample_id = c("s1", "s2"),
      stringsAsFactors = FALSE
    ),
    updateSelectInputFn = function(session, inputId, choices, selected) {
      design_calls[[length(design_calls) + 1]] <<- list(
        session = session,
        inputId = inputId,
        choices = choices,
        selected = selected
      )
    },
    messageFn = function(text) {
      design_messages <<- c(design_messages, text)
    }
  )
  expect_length(design_calls, 3)
  expect_equal(design_calls[[1]]$inputId, "color_variable")
  expect_equal(design_calls[[3]]$inputId, "ruv_grouping_variable")
  expect_true(any(grepl("Updated grouping variable choices", design_messages, fixed = TRUE)))

  label_plot <- build_label_plot("a)")
  title_plot <- build_title_plot("Pre-Normalisation")
  expect_s3_class(label_plot, "ggplot")
  expect_s3_class(title_plot, "ggplot")

  img_path <- tempfile(fileext = ".png")
  png::writePNG(array(1, dim = c(1, 1, 4)), img_path)
  image_plot <- load_image_plot(img_path, messageFn = function(...) NULL)
  missing_plot <- load_image_plot(file.path(tempdir(), "missing-prot-norm-image.png"), messageFn = function(...) NULL)
  expect_s3_class(image_plot, "ggplot")
  expect_s3_class(missing_plot, "ggplot")
})

test_that("prot norm shared normalization and RUV helpers preserve workflow building blocks", {
  skipIfMissingMultiScholaRBindings(
    "prepareProtNormNormalizationRun",
    "runProtNormBetweenSamplesStep",
    "runProtNormPostNormalizationQcStep",
    "buildProtNormSkippedRuvResult",
    "applyProtNormSkippedRuvState",
    "buildProtNormManualRuvResult",
    "updateProtNormRuvAuditTrail",
    "persistProtNormRuvResult"
  )

  prepare_run <- getMultiScholaRBinding("prepareProtNormNormalizationRun")
  run_between_samples <- getMultiScholaRBinding("runProtNormBetweenSamplesStep")
  run_post_qc <- getMultiScholaRBinding("runProtNormPostNormalizationQcStep")
  build_skipped_result <- getMultiScholaRBinding("buildProtNormSkippedRuvResult")
  apply_skipped_state <- getMultiScholaRBinding("applyProtNormSkippedRuvState")
  build_manual_result <- getMultiScholaRBinding("buildProtNormManualRuvResult")
  update_audit_trail <- getMultiScholaRBinding("updateProtNormRuvAuditTrail")
  persist_ruv_result <- getMultiScholaRBinding("persistProtNormRuvResult")

  current_s4 <- makeProtNormCharacterizationData(label = "norm-current")
  norm_data <- new.env(parent = emptyenv())
  run_context <- prepare_run(
    stateManager = list(
      current_state = "protein_replicate_filtered",
      getState = function(state) current_s4
    ),
    normData = norm_data,
    reqFn = function(x) x,
    initialiseGridFn = function() "grid_obj",
    messageFn = function(...) NULL
  )
  expect_equal(run_context$currentState, "protein_replicate_filtered")
  expect_identical(run_context$currentS4, current_s4)
  expect_equal(norm_data$QC_composite_figure, "grid_obj")

  saved_matrix <- NULL
  checkpoint_call <- NULL
  normalized_s4 <- run_between_samples(
    currentS4 = current_s4,
    normMethod = "cyclicloess",
    normData = norm_data,
    proteinQcDir = tempdir(),
    normaliseBetweenSamplesFn = function(currentS4, normalisation_method) {
      expect_equal(normalisation_method, "cyclicloess")
      currentS4
    },
    captureCheckpointFn = function(object, checkpoint, label) {
      checkpoint_call <<- list(object = object, checkpoint = checkpoint, label = label)
    },
    existsFn = function(...) FALSE,
    saveMatrixFn = function(data, path) {
      saved_matrix <<- list(data = data, path = path)
    },
    messageFn = function(...) NULL
  )
  expect_identical(normalized_s4, current_s4)
  expect_identical(norm_data$normalized_protein_obj, current_s4)
  expect_equal(checkpoint_call$checkpoint, "cp05")
  expect_true(grepl("normalized_protein_matrix_pre_ruv.tsv$", saved_matrix$path))

  post_qc_calls <- list()
  run_post_qc(
    normalizedS4 = current_s4,
    normData = norm_data,
    generatePostNormalizationQcFn = function(object) {
      post_qc_calls[[length(post_qc_calls) + 1]] <<- object
    },
    messageFn = function(...) NULL
  )
  expect_true(isTRUE(norm_data$normalization_complete))
  expect_identical(post_qc_calls[[1]], current_s4)

  skipped_result <- build_skipped_result()
  expect_true(isTRUE(skipped_result$ruv_skipped))
  expect_match(skipped_result$skip_reason, "User selected skip", fixed = TRUE)

  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- new.env(parent = emptyenv())
  saved_state <- NULL
  workflow_data$state_manager$saveState <- function(...) {
    saved_state <<- list(...)
  }
  saved_rds <- NULL
  apply_skipped_state(
    normalizedS4 = current_s4,
    normMethod = "cyclicloess",
    normData = norm_data,
    workflowData = workflow_data,
    sourceDir = tempdir(),
    skipResult = skipped_result,
    saveRdsFn = function(object, path) {
      saved_rds <<- list(object = object, path = path)
    },
    messageFn = function(...) NULL
  )
  expect_identical(norm_data$ruv_normalized_obj, current_s4)
  expect_true(isTRUE(norm_data$ruv_complete))
  expect_true(is.na(norm_data$best_k))
  expect_equal(saved_state$state_name, "normalized")
  expect_true(grepl("ruv_optimization_results.RDS$", saved_rds$path))

  manual_result <- build_manual_result(
    percentageAsNegCtrl = 5,
    ruvK = 2,
    controlGenesIndex = c(TRUE, FALSE, TRUE),
    cancorPlot = "plot_obj"
  )
  expect_equal(manual_result$best_k, 2)
  expect_equal(manual_result$optimization_results$num_controls, 2)

  updated_config <- NULL
  update_audit_trail(
    ruvK = 3,
    controlGenesIndex = c(TRUE, FALSE),
    percentageAsNegCtrl = 10,
    modeLabel = "automatic",
    existsFn = function(...) TRUE,
    getFn = function(...) list(existing = TRUE),
    assignFn = function(name, value, envir) {
      updated_config <<- list(name = name, value = value)
    },
    updateRuvParametersFn = function(config_list, ruv_k, control_genes_index, percentage_as_neg_ctrl) {
      config_list$ruv_k <- ruv_k
      config_list$num_controls <- sum(control_genes_index, na.rm = TRUE)
      config_list$percentage <- percentage_as_neg_ctrl
      config_list
    },
    messageFn = function(...) NULL
  )
  expect_equal(updated_config$name, "config_list")
  expect_equal(updated_config$value$ruv_k, 3)
  expect_equal(updated_config$value$num_controls, 1)

  persisted_messages <- character()
  persisted_result <- persist_ruv_result(
    ruvResult = manual_result,
    workflowData = workflow_data,
    sourceDir = tempdir(),
    resultLabel = "manual RUV results",
    saveRdsFn = function(object, path) {
      persisted_messages <<- c(persisted_messages, path)
    },
    catFn = function(...) NULL
  )
  expect_identical(persisted_result, manual_result)
  expect_identical(workflow_data$ruv_optimization_result, manual_result)
  expect_true(any(grepl("ruv_optimization_results.RDS$", persisted_messages)))
})

test_that("prot norm shared RUV cleanup and composite helpers preserve orchestration", {
  skipIfMissingMultiScholaRBindings(
    "resolveProtNormRuvParameters",
    "applyProtNormRuvCorrectionStep",
    "finalizeProtNormRuvCleanupStep",
    "resolveProtNormStep6QcObject",
    "runProtNormStep6RuvQc",
    "resolveProtNormCompositeFigureInputs",
    "generateProtNormCompositeQcFigure",
    "finalizeProtNormWorkflowState",
    "buildProtNormCompletionNotification"
  )

  resolve_ruv_params <- getMultiScholaRBinding("resolveProtNormRuvParameters")
  apply_ruv_step <- getMultiScholaRBinding("applyProtNormRuvCorrectionStep")
  finalize_cleanup <- getMultiScholaRBinding("finalizeProtNormRuvCleanupStep")
  resolve_step6_qc <- getMultiScholaRBinding("resolveProtNormStep6QcObject")
  run_step6_qc <- getMultiScholaRBinding("runProtNormStep6RuvQc")
  resolve_composite_inputs <- getMultiScholaRBinding("resolveProtNormCompositeFigureInputs")
  generate_composite_qc <- getMultiScholaRBinding("generateProtNormCompositeQcFigure")
  finalize_workflow <- getMultiScholaRBinding("finalizeProtNormWorkflowState")
  build_completion_notification <- getMultiScholaRBinding("buildProtNormCompletionNotification")

  normalized_s4 <- makeProtNormCharacterizationData(label = "ruv-source")
  norm_data <- new.env(parent = emptyenv())
  workflow_data <- new.env(parent = emptyenv())

  automatic_result <- resolve_ruv_params(
    normalizedS4 = normalized_s4,
    input = list(
      ruv_mode = "automatic",
      auto_percentage_min = 1,
      auto_percentage_max = 3,
      separation_metric = "max_difference",
      k_penalty_weight = 0.5,
      adaptive_k_penalty = TRUE
    ),
    normData = norm_data,
    workflowData = workflow_data,
    sourceDir = tempdir(),
    getRuvGroupingVariableFn = function() "group",
    withProgressFn = function(message, detail, value, expr) force(expr),
    findBestNegCtrlPercentageFn = function(...) list(
      best_percentage = 2,
      best_k = 3,
      best_control_genes_index = c(TRUE, FALSE),
      best_composite_score = 0.8,
      best_separation_score = 0.6,
      optimization_results = data.frame(percentage = 1:3)
    ),
    persistRuvResultFn = function(ruvResult, workflowData, sourceDir, resultLabel, ...) {
      workflowData$ruv_optimization_result <- ruvResult
      invisible(ruvResult)
    },
    updateAuditTrailFn = function(...) NULL,
    messageFn = function(...) NULL
  )
  expect_equal(automatic_result$percentageAsNegCtrl, 2)
  expect_equal(norm_data$best_k, 3)
  expect_equal(sum(norm_data$control_genes_index, na.rm = TRUE), 1)

  corrected_s4 <- apply_ruv_step(
    normalizedS4 = normalized_s4,
    normData = norm_data,
    getRuvGroupingVariableFn = function() "group",
    ruvIII_C_VaryingFn = function(object, ruv_grouping_variable, ruv_number_k, ctrl) {
      expect_identical(object, normalized_s4)
      expect_equal(ruv_grouping_variable, "group")
      expect_equal(ruv_number_k, 3)
      normalized_s4
    },
    captureCheckpointFn = function(...) NULL,
    messageFn = function(...) NULL
  )
  expect_identical(corrected_s4, normalized_s4)
  expect_identical(norm_data$ruv_normalized_obj, normalized_s4)

  workflow_data$state_manager <- new.env(parent = emptyenv())
  saved_state <- NULL
  workflow_data$state_manager$saveState <- function(...) saved_state <<- list(...)
  workflow_data$protein_counts <- list()
  cleaned_s4 <- finalize_cleanup(
    ruvCorrectedS4 = normalized_s4,
    input = list(norm_method = "cyclicloess", ruv_mode = "automatic", ruv_percentage = 5),
    normData = norm_data,
    workflowData = workflow_data,
    omicType = "proteomics",
    experimentLabel = "exp1",
    removeRowsWithMissingValuesPercentFn = function(theObject) theObject,
    countDistinctProteinsFn = function(s4Object) 2,
    updateProteinFilteringFn = function(...) "filter_plot",
    messageFn = function(...) NULL
  )
  expect_identical(cleaned_s4, normalized_s4)
  expect_equal(workflow_data$protein_counts$after_ruv_filtering, 2)
  expect_equal(norm_data$post_norm_filtering_plot, "filter_plot")
  expect_equal(saved_state$state_name, "ruv_corrected")

  fallback_norm_data <- new.env(parent = emptyenv())
  fallback_norm_data$ruv_normalized_obj <- normalized_s4
  expect_identical(
    resolve_step6_qc(step5Object = NULL, normData = fallback_norm_data, messageFn = function(...) NULL),
    normalized_s4
  )

  fallback_norm_data$ruv_optimization_result <- list(best_cancor_plot = "plot_obj")
  run_step6_qc(
    ruvMode = "automatic",
    step6Object = normalized_s4,
    normData = fallback_norm_data,
    qcDir = tempdir(),
    generateRuvCorrectedQcFn = function(object) expect_identical(object, normalized_s4),
    ggsaveFn = function(...) NULL,
    initPathsFn = function(paths) getMultiScholaRBinding("initializeProtNormQcPlotPaths")(paths),
    messageFn = function(...) NULL
  )
  expect_true(grepl("ruv_corrected_cancor.png$", fallback_norm_data$qc_plot_paths$ruv_corrected$cancor))

  skip_inputs <- resolve_composite_inputs("skip", tempdir())
  full_inputs <- resolve_composite_inputs("automatic", tempdir())
  expect_equal(skip_inputs$ncol, 2)
  expect_equal(full_inputs$ncol, 3)
  expect_length(full_inputs$plotFiles, 15)

  composite_saved <- list()
  generate_composite_qc(
    ruvMode = "automatic",
    qcDir = tempdir(),
    omicType = "proteomics",
    resolveInputsFn = function(...) list(
      ncol = 3,
      plotFiles = "plot.png",
      rowLabels = list(pca = "a)"),
      columnLabels = c("Pre")
    ),
    generateCompositeFn = function(...) list(plot = "composite_plot", width = 12, height = 8),
    savePlotFn = function(plot, base_path, plot_name, width, height, ...) {
      composite_saved <<- list(plot = plot, base_path = base_path, plot_name = plot_name, width = width, height = height)
    },
    messageFn = function(...) NULL
  )
  expect_equal(composite_saved$plot_name, "proteomics_composite_QC_figure")

  finalize_workflow(
    normData = fallback_norm_data,
    gcFn = function() NULL,
    messageFn = function(...) NULL
  )
  expect_true(isTRUE(fallback_norm_data$ruv_complete))
  expect_match(build_completion_notification("skip"), "RUV skipped", fixed = TRUE)
  expect_match(build_completion_notification("automatic"), "Normalization and RUV correction completed", fixed = TRUE)
})

test_that("prot norm shared export and reset helpers preserve session handoff", {
  skipIfMissingMultiScholaRBindings(
    "canProtNormExportFilteredSession",
    "resolveProtNormExportSourceDir",
    "collectProtNormExportSessionData",
    "resolveProtNormPreNormalizationState",
    "revertProtNormStateManagerToPreNormalization",
    "resetProtNormOutputs"
  )

  can_export_session <- getMultiScholaRBinding("canProtNormExportFilteredSession")
  resolve_source_dir <- getMultiScholaRBinding("resolveProtNormExportSourceDir")
  collect_session_data <- getMultiScholaRBinding("collectProtNormExportSessionData")
  resolve_pre_state <- getMultiScholaRBinding("resolveProtNormPreNormalizationState")
  revert_pre_state <- getMultiScholaRBinding("revertProtNormStateManagerToPreNormalization")
  reset_outputs <- getMultiScholaRBinding("resetProtNormOutputs")

  expect_true(can_export_session(TRUE, "obj"))
  expect_false(can_export_session(FALSE, "obj"))
  expect_error(
    resolve_source_dir(list(source_dir = file.path(tempdir(), "missing")), dirExistsFn = function(path) FALSE),
    "Could not find the source directory"
  )
  expect_equal(
    resolve_source_dir(list(source_dir = tempdir()), dirExistsFn = function(path) TRUE),
    tempdir()
  )

  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- new.env(parent = emptyenv())
  workflow_data$state_manager$current_state <- "correlation_filtered"
  workflow_data$state_manager$getState <- function(state) makeProtNormCharacterizationData(label = state)
  workflow_data$contrasts_tbl <- data.frame(friendly_names = "A vs B", stringsAsFactors = FALSE)
  workflow_data$design_matrix <- data.frame(group = c("A", "B"), stringsAsFactors = FALSE)
  workflow_data$config_list <- list(globalParameters = list(workflow_type = "DIA"))
  workflow_data$fasta_metadata <- list(db = "ref")
  workflow_data$accession_cleanup_results <- list(clean = TRUE)
  workflow_data$ruv_optimization_result <- list(best_k = 2)
  workflow_data$qc_params <- list(minimum = 1)
  workflow_data$protein_counts <- list(final_for_de = 2)
  workflow_data$mixed_species_analysis <- list(enabled = TRUE)
  norm_data <- new.env(parent = emptyenv())
  norm_data$correlation_filtered_obj <- makeProtNormCharacterizationData(label = "correlation_filtered")
  norm_data$best_k <- 2
  norm_data$correlation_threshold <- 0.75

  session_data <- collect_session_data(
    workflowData = workflow_data,
    normData = norm_data,
    input = list(norm_method = "cyclicloess", ruv_mode = "automatic"),
    timeFn = function() as.POSIXct("2026-04-11 10:00:00", tz = "UTC"),
    messageFn = function(...) NULL
  )
  expect_equal(session_data$normalization_method, "cyclicloess")
  expect_equal(session_data$ruv_mode, "automatic")
  expect_equal(session_data$final_protein_count, 2)
  expect_equal(session_data$final_sample_count, 2)

  expect_equal(
    resolve_pre_state(c("raw_data_s4", "sample_filtered", "protein_replicate_filtered", "normalized")),
    "protein_replicate_filtered"
  )

  workflow_data$state_manager$getHistory <- function() c("raw_data_s4", "sample_filtered", "protein_replicate_filtered", "normalized")
  workflow_data$state_manager$revertToState <- function(state) makeProtNormCharacterizationData(label = state)
  revert_result <- revert_pre_state(
    workflowData = workflow_data,
    messageFn = function(...) NULL
  )
  expect_equal(revert_result$previousState, "protein_replicate_filtered")
  expect_s4_class(revert_result$revertedS4, "mockProtNormCharacterizationData")

  output <- new.env(parent = emptyenv())
  reset_outputs(
    output = output,
    ruvMode = "automatic",
    groupingVariable = "group",
    renderTextFn = function(expr) force(expr)
  )
  expect_match(output$correlation_filter_summary, "No correlation filtering applied yet", fixed = TRUE)
  expect_match(output$filtering_summary_text, "Filtering summary will be available", fixed = TRUE)
  expect_match(output$ruv_optimization_summary, "Run normalization to see RUV optimization results", fixed = TRUE)
})
