module_ci_prot_norm_design <- function(samples = c("A1", "A2", "B1", "B2"),
                                       groups = c("A", "A", "B", "B"),
                                       replicates = c("R1", "R1", "R2", "R2"),
                                       batches = c("Batch1", "Batch2", "Batch1", "Batch2")) {
  groups <- rep(groups, length.out = length(samples))
  replicates <- rep(replicates, length.out = length(samples))
  batches <- rep(batches, length.out = length(samples))

  data.frame(
    Run = samples,
    group = groups,
    replicates = replicates,
    batch = batches,
    stringsAsFactors = FALSE
  )
}

module_ci_prot_norm_values <- function(samples = c("A1", "A2", "B1", "B2")) {
  values <- rbind(
    P_STABLE = c(10, 11, 12, 13),
    P_SHIFT = c(100, 110, 120, 130),
    P_ZERO = c(0, 0, 0, 0),
    P_NEG = c(-1, -2, -3, -4),
    P_NA = c(NA, 5, NA, 7),
    P_INF = c(Inf, 1, 2, 3),
    P_NAN = c(NaN, 4, 5, 6),
    P_CONST = c(42, 42, 42, 42),
    P_CTRL_1 = c(20.0, 20.1, 20.2, 20.1),
    P_CTRL_2 = c(21.0, 21.1, 21.2, 21.1),
    P_CTRL_3 = c(22.0, 22.1, 22.2, 22.1),
    P_CTRL_4 = c(23.0, 23.1, 23.2, 23.1),
    P_CTRL_5 = c(24.0, 24.1, 24.2, 24.1),
    P_CTRL_6 = c(25.0, 25.1, 25.2, 25.1)
  )
  colnames(values) <- samples
  values
}

module_ci_prot_norm_table <- function(values = module_ci_prot_norm_values()) {
  output <- data.frame(
    Protein.Ids = rownames(values),
    as.data.frame(values, check.names = FALSE),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  rownames(output) <- NULL
  output
}

module_ci_prot_norm_control_index <- function(features = rownames(module_ci_prot_norm_values()),
                                              selected = grep("^P_CTRL_", features, value = TRUE)) {
  controls <- features %in% selected
  names(controls) <- features
  controls
}

module_ci_prot_norm_args <- function(features = rownames(module_ci_prot_norm_values()),
                                     workflow_type = "DIA",
                                     report_template = NULL,
                                     use_limpa = FALSE,
                                     normalisation_method = "none",
                                     ruv_grouping_variable = "group",
                                     ruv_number_k = 2L,
                                     ctrl = module_ci_prot_norm_control_index(features)) {
  list(
    globalParameters = list(
      workflow_type = workflow_type,
      report_template = report_template,
      use_limpa = use_limpa
    ),
    normaliseBetweenSamples = list(
      normalisation_method = normalisation_method
    ),
    getNegCtrlProtAnova = list(
      ruv_grouping_variable = ruv_grouping_variable,
      percentage_as_neg_ctrl = 40,
      num_neg_ctrl = sum(ctrl, na.rm = TRUE),
      ruv_qval_cutoff = 0.05,
      ruv_fdr_method = "BH",
      exclude_pool_samples = TRUE
    ),
    ruvCancor = list(
      ctrl = ctrl,
      num_components_to_impute = 2L,
      ruv_grouping_variable = ruv_grouping_variable
    ),
    ruvIII_C_Varying = list(
      ruv_grouping_variable = ruv_grouping_variable,
      ruv_number_k = ruv_number_k,
      ctrl = ctrl
    ),
    removeRowsWithMissingValuesPercent = list(
      ruv_grouping_variable = ruv_grouping_variable,
      groupwise_percentage_cutoff = 50,
      max_groups_percentage_cutoff = 50,
      proteins_intensity_cutoff_percentile = 1
    ),
    filterSamplesByProteinCorrelationThreshold = list(
      min_pearson_correlation_threshold = 0.75
    )
  )
}

module_ci_prot_norm_object <- function(protein_quant_table = module_ci_prot_norm_table(),
                                       design = module_ci_prot_norm_design(colnames(protein_quant_table)[-1]),
                                       args = module_ci_prot_norm_args(protein_quant_table$Protein.Ids)) {
  ProteinQuantitativeData(
    protein_quant_table = protein_quant_table,
    protein_id_column = "Protein.Ids",
    protein_id_table = dplyr::distinct(protein_quant_table, Protein.Ids),
    design_matrix = design,
    sample_id = "Run",
    group_id = "group",
    technical_replicate_id = "replicates",
    args = args
  )
}

module_ci_prot_norm_matrix <- function(object) {
  mat <- as.matrix(object@protein_quant_table[, setdiff(colnames(object@protein_quant_table), object@protein_id_column), drop = FALSE])
  rownames(mat) <- object@protein_quant_table[[object@protein_id_column]]
  mat
}

module_ci_prot_norm_assert_sanity <- function(object, allow_na = TRUE) {
  mat <- module_ci_prot_norm_matrix(object)
  testthat::expect_false(any(is.nan(mat)))
  testthat::expect_false(any(is.infinite(mat)))
  if (!allow_na) {
    testthat::expect_false(any(is.na(mat)))
  }
  testthat::expect_setequal(
    setdiff(colnames(object@protein_quant_table), object@protein_id_column),
    as.character(object@design_matrix[[object@sample_id]])
  )
  invisible(TRUE)
}

module_ci_prot_norm_make_function_with_overrides <- function(fun, replacements) {
  fun_override <- fun
  environment(fun_override) <- list2env(replacements, parent = environment(fun))
  fun_override
}

module_ci_prot_norm_param_reader <- function(theObject, name, default = NULL) {
  param_value <- tryCatch(get(name, envir = parent.frame(), inherits = FALSE), error = function(e) NULL)
  if (!is.null(param_value)) {
    return(param_value)
  }

  for (section_name in names(theObject@args)) {
    section <- theObject@args[[section_name]]
    if (is.list(section) && name %in% names(section) && !is.null(section[[name]])) {
      return(section[[name]])
    }
  }

  default
}

module_ci_prot_norm_param_updater <- function(theObject, name) {
  param_value <- tryCatch(get(name, envir = parent.frame(), inherits = FALSE), error = function(e) NULL)
  if (is.null(theObject@args$module_ci_param_audit)) {
    theObject@args$module_ci_param_audit <- list()
  }
  theObject@args$module_ci_param_audit[[name]] <- param_value
  theObject
}

module_ci_prot_norm_s4_methods <- function() {
  normalizers <- list(
    normalizeCyclicLoess = function(x) x + 10,
    normalizeQuantiles = function(x) x + 20,
    normalizeMedianAbsValues = function(x) x + 30
  )

  overrides <- list(
    checkParamsObjectFunctionSimplify = module_ci_prot_norm_param_reader,
    updateParamInObject = module_ci_prot_norm_param_updater,
    cleanDesignMatrix = function(theObject) theObject,
    resolveProtLimmaNormalizer = function(functionName, envir = parent.frame()) normalizers[[functionName]],
    impute.nipals = function(x, ncomp) {
      x[is.na(x)] <- -99
      x
    },
    ruv_cancorplot = function(Y, X, ctl) {
      list(Y = Y, X = X, ctl = ctl)
    },
    getRuvIIIReplicateMatrixHelper = function(design_matrix, sample_id_column, grouping_variable) {
      samples <- as.character(design_matrix$Run)
      replicate_matrix <- diag(length(samples), nrow = length(samples), ncol = length(samples))
      dimnames(replicate_matrix) <- list(samples, samples)
      replicate_matrix
    },
    RUVIII_C_Varying = function(k, Y, M, toCorrect, potentialControls) {
      corrected <- Y + as.numeric(k) / 10
      dimnames(corrected) <- dimnames(Y)
      corrected
    }
  )

  list(
    normalise = module_ci_prot_norm_make_function_with_overrides(
      methods::selectMethod("normaliseBetweenSamples", "ProteinQuantitativeData"),
      overrides
    ),
    cancor = module_ci_prot_norm_make_function_with_overrides(
      methods::selectMethod("ruvCancor", "ProteinQuantitativeData"),
      overrides
    ),
    ruviii = module_ci_prot_norm_make_function_with_overrides(
      methods::selectMethod("ruvIII_C_Varying", "ProteinQuantitativeData"),
      overrides
    )
  )
}

module_ci_prot_norm_cancor_plot <- function(deltas = c(0.05, 0.20, 0.19, 0.18)) {
  plot_data <- data.frame(
    featureset = rep(c("Control", "All"), each = length(deltas)),
    K = rep(seq_along(deltas), times = 2),
    cc = c(rep(0.10, length(deltas)), 0.10 + deltas)
  )
  ggplot2::ggplot(plot_data, ggplot2::aes(K, cc, colour = featureset)) +
    ggplot2::geom_point()
}

module_ci_prot_norm_state_manager <- function(states,
                                              current_state = names(states)[[length(states)]],
                                              history = names(states)) {
  manager <- new.env(parent = emptyenv())
  manager$states <- states
  manager$state_history <- history
  manager$current_state <- current_state
  manager$getState <- function(state_name = NULL) {
    if (is.null(state_name)) {
      state_name <- manager$current_state
    }
    manager$states[[state_name]]
  }
  manager$saveState <- function(state_name, s4_data_object, config_object = list(), description = "") {
    manager$states[[state_name]] <- s4_data_object
    manager$state_history <- c(manager$state_history, state_name)
    manager$current_state <- state_name
    invisible(TRUE)
  }
  manager$getHistory <- function() manager$state_history
  manager$revertToState <- function(state_name) {
    if (is.null(manager$states[[state_name]])) {
      stop(sprintf("State '%s' not found.", state_name), call. = FALSE)
    }
    manager$current_state <- state_name
    manager$states[[state_name]]
  }
  manager
}

module_ci_prot_norm_workflow <- function(current_object = module_ci_prot_norm_object(),
                                         workflow_type = "DIA",
                                         report_template = NULL,
                                         use_limpa = FALSE,
                                         ruv_result = list(best_k = 2L, best_percentage = 40)) {
  workflow <- new.env(parent = emptyenv())
  workflow$state_manager <- module_ci_prot_norm_state_manager(
    states = list(
      protein_replicate_filtered = current_object,
      normalized = current_object,
      ruv_corrected = current_object,
      correlation_filtered = current_object
    ),
    current_state = "correlation_filtered",
    history = c("protein_replicate_filtered", "normalized", "ruv_corrected", "correlation_filtered")
  )
  workflow$design_matrix <- current_object@design_matrix
  workflow$contrasts_tbl <- data.frame(
    contrasts = "groupA-groupB",
    full_format = "A_vs_B=groupA-groupB",
    friendly_names = "A_vs_B",
    stringsAsFactors = FALSE
  )
  workflow$config_list <- list(
    globalParameters = list(
      workflow_type = workflow_type,
      report_template = report_template,
      use_limpa = use_limpa
    )
  )
  workflow$qc_params <- list(normalization = list(methods = c("none", "cyclicloess", "quantile", "scale")))
  workflow$protein_counts <- list(after_qc_filtering = nrow(current_object@protein_quant_table))
  workflow$ruv_optimization_result <- ruv_result
  workflow$tab_status <- list(normalization = "pending", differential_expression = "disabled")
  workflow$state_update_trigger <- NULL
  workflow
}

module_ci_prot_norm_data <- function(current_object = module_ci_prot_norm_object(),
                                     ruv_result = list(best_k = 2L, best_percentage = 40)) {
  norm_data <- new.env(parent = emptyenv())
  norm_data$pre_norm_qc_generated <- FALSE
  norm_data$normalization_complete <- FALSE
  norm_data$ruv_complete <- FALSE
  norm_data$correlation_filtering_complete <- FALSE
  norm_data$qc_plot_paths <- NULL
  norm_data$normalized_protein_obj <- current_object
  norm_data$ruv_normalized_obj <- current_object
  norm_data$correlation_filtered_obj <- current_object
  norm_data$best_k <- ruv_result$best_k
  norm_data$control_genes_index <- module_ci_prot_norm_control_index(current_object@protein_quant_table$Protein.Ids)
  norm_data$correlation_vector <- NULL
  norm_data$correlation_threshold <- 0.75
  norm_data$final_qc_plot <- NULL
  norm_data$final_filtering_plot <- NULL
  norm_data$post_norm_filtering_plot <- NULL
  norm_data$filtering_summary_text <- NULL
  norm_data$ruv_optimization_result <- ruv_result
  norm_data
}

module_ci_prot_norm_input <- function(norm_method = "none",
                                      ruv_mode = "manual",
                                      ruv_percentage = 40,
                                      ruv_k = 2L,
                                      min_threshold = 0.75) {
  list(
    norm_method = norm_method,
    ruv_mode = ruv_mode,
    ruv_percentage = ruv_percentage,
    ruv_k = ruv_k,
    auto_percentage_min = 10,
    auto_percentage_max = 12,
    separation_metric = "max_difference",
    k_penalty_weight = 0.5,
    adaptive_k_penalty = TRUE,
    min_pearson_correlation_threshold = min_threshold
  )
}

module_ci_prot_norm_correlation_pairs <- function(case = c("pass_all", "fail_one", "fail_many", "boundary", "missing_sample")) {
  case <- match.arg(case)
  switch(case,
    pass_all = data.frame(
      Run.x = c("A2", "B2"),
      Run.y = c("A1", "B1"),
      pearson_correlation = c(0.99, 0.98),
      stringsAsFactors = FALSE
    ),
    fail_one = data.frame(
      Run.x = c("A2", "B2"),
      Run.y = c("A1", "B1"),
      pearson_correlation = c(0.99, 0.20),
      stringsAsFactors = FALSE
    ),
    fail_many = data.frame(
      Run.x = c("A2", "B2"),
      Run.y = c("A1", "B1"),
      pearson_correlation = c(0.20, 0.30),
      stringsAsFactors = FALSE
    ),
    boundary = data.frame(
      Run.x = c("A2", "B2"),
      Run.y = c("A1", "B1"),
      pearson_correlation = c(0.75, 0.75),
      stringsAsFactors = FALSE
    ),
    missing_sample = data.frame(
      Run.x = "Ghost",
      Run.y = "A1",
      pearson_correlation = 0.99,
      stringsAsFactors = FALSE
    )
  )
}

module_ci_prot_norm_empty_pairs <- function() {
  data.frame(
    Run.x = character(),
    Run.y = character(),
    pearson_correlation = numeric(),
    stringsAsFactors = FALSE
  )
}

module_ci_prot_norm_expected_qc_names <- function(prefix) {
  paste0(prefix, c("_pca.png", "_rle.png", "_density.png", "_correlation.png"))
}

module_ci_prot_norm_artifact_saver <- function(qcDir, filename, plotObject, width, height, dpi = 150) {
  dir.create(qcDir, recursive = TRUE, showWarnings = FALSE)
  artifact_path <- file.path(qcDir, filename)
  writeLines("module-ci plot placeholder", artifact_path)
  artifact_path
}

module_ci_prot_norm_ggplot_stub <- function(label = "plot") {
  ggplot2::ggplot(data.frame(x = 1, y = 1, label = label), ggplot2::aes(x, y)) +
    ggplot2::geom_point() +
    ggplot2::ggtitle(label)
}

module_ci_prot_norm_export_variant <- function(ruv_mode = "manual",
                                               workflow_type = "DIA",
                                               report_template = NULL,
                                               use_limpa = FALSE,
                                               norm_method = "quantile") {
  ruv_result <- if (identical(ruv_mode, "skip")) {
    buildProtNormSkippedRuvResult()
  } else {
    list(
      best_k = 2L,
      best_percentage = 40,
      best_control_genes_index = module_ci_prot_norm_control_index(),
      optimization_results = data.frame(
        percentage_requested = 40,
        best_k = 2L,
        composite_score = 0.5,
        status = "ok",
        stringsAsFactors = FALSE
      )
    )
  }

  object <- module_ci_prot_norm_object(args = module_ci_prot_norm_args(
    workflow_type = workflow_type,
    report_template = report_template,
    use_limpa = use_limpa
  ))
  workflow <- module_ci_prot_norm_workflow(
    current_object = object,
    workflow_type = workflow_type,
    report_template = report_template,
    use_limpa = use_limpa,
    ruv_result = ruv_result
  )
  norm_data <- module_ci_prot_norm_data(object, ruv_result = ruv_result)
  input <- module_ci_prot_norm_input(norm_method = norm_method, ruv_mode = ruv_mode)
  collectProtNormExportSessionData(
    workflowData = workflow,
    normData = norm_data,
    input = input,
    timeFn = function() as.POSIXct("2026-05-05 12:00:00", tz = "UTC"),
    messageFn = function(...) invisible(NULL)
  )
}
