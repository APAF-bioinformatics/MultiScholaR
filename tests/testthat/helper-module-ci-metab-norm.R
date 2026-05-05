module_ci_metab_norm_samples <- function() {
  c("S1", "S2", "S3", "S4")
}

module_ci_metab_norm_design <- function(samples = module_ci_metab_norm_samples(),
                                        small_n = FALSE) {
  groups <- if (isTRUE(small_n)) {
    c("WT", "KO", "RES", "REC")[seq_along(samples)]
  } else {
    rep(c("WT", "KO"), each = 2, length.out = length(samples))
  }

  data.frame(
    Run = samples,
    group = groups,
    batch = rep(c("B1", "B2"), length.out = length(samples)),
    tech_rep_group = rep(c("Pair1", "Pair1", "Pair2", "Pair2"), length.out = length(samples)),
    stringsAsFactors = FALSE
  )
}

module_ci_metab_norm_values <- function(samples = module_ci_metab_norm_samples(),
                                        positive_only = FALSE,
                                        assay_shift = 0) {
  values <- rbind(
    M_STABLE = c(100, 110, 120, 130),
    M_ZERO = c(0, 0, 50, 60),
    M_MISSING = c(NA, 80, 85, 90),
    M_NEG = c(-5, 10, 15, 20),
    M_CONST = c(42, 42, 42, 42),
    M_SHIFT = c(200, 210, 230, 260),
    CTRL_1 = c(20.0, 20.2, 20.4, 20.6),
    CTRL_2 = c(21.0, 21.2, 21.4, 21.6),
    CTRL_3 = c(22.0, 22.2, 22.4, 22.6),
    CTRL_4 = c(23.0, 23.2, 23.4, 23.6),
    CTRL_5 = c(24.0, 24.2, 24.4, 24.6),
    CTRL_6 = c(25.0, 25.2, 25.4, 25.6),
    IS_A = c(1000, 900, 1100, 1000),
    IS_B = c(2000, 2100, 1900, 2050)
  )

  values <- values[, seq_along(samples), drop = FALSE]
  colnames(values) <- samples
  values[!is.na(values)] <- values[!is.na(values)] + assay_shift

  if (isTRUE(positive_only)) {
    values[is.na(values)] <- 5
    values[values <= 0] <- abs(values[values <= 0]) + 1
  }

  values
}

module_ci_metab_norm_assay <- function(samples = module_ci_metab_norm_samples(),
                                       include_itsd = TRUE,
                                       include_controls = TRUE,
                                       positive_only = FALSE,
                                       assay_shift = 0) {
  values <- module_ci_metab_norm_values(
    samples = samples,
    positive_only = positive_only,
    assay_shift = assay_shift
  )

  keep <- rep(TRUE, nrow(values))
  if (!isTRUE(include_itsd)) {
    keep <- keep & !grepl("^IS_", rownames(values))
  }
  if (!isTRUE(include_controls)) {
    keep <- keep & !grepl("^CTRL_", rownames(values))
  }
  values <- values[keep, , drop = FALSE]

  sample_df <- as.data.frame(values, check.names = FALSE)
  names(sample_df) <- samples

  data.frame(
    database_identifier = rownames(values),
    annotation_id = rownames(values),
    Name = rownames(values),
    sample_df,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

module_ci_metab_norm_assays <- function(layout = "combined",
                                        samples = module_ci_metab_norm_samples(),
                                        include_itsd = TRUE,
                                        include_controls = TRUE,
                                        positive_only = FALSE) {
  primary <- module_ci_metab_norm_assay(
    samples = samples,
    include_itsd = include_itsd,
    include_controls = include_controls,
    positive_only = positive_only,
    assay_shift = 0
  )
  secondary <- module_ci_metab_norm_assay(
    samples = samples,
    include_itsd = include_itsd,
    include_controls = include_controls,
    positive_only = positive_only,
    assay_shift = 5
  )

  switch(
    layout,
    lc = list(LCMS_Pos = primary),
    gc = list(GCMS = primary),
    combined = list(LCMS_Pos = primary, GCMS = secondary),
    list(LCMS_Pos = primary)
  )
}

module_ci_metab_norm_args <- function(normalisation_method = "none",
                                      formula_string = "~ 0 + group",
                                      ruv_grouping_variable = "tech_rep_group") {
  list(
    normaliseBetweenSamples = list(normalisation_method = normalisation_method),
    daAnalysisParameters = list(formula_string = formula_string),
    ruvIII_C_Varying = list(ruv_grouping_variable = ruv_grouping_variable)
  )
}

module_ci_metab_norm_object <- function(layout = "combined",
                                        samples = module_ci_metab_norm_samples(),
                                        include_itsd = TRUE,
                                        include_controls = TRUE,
                                        positive_only = FALSE,
                                        args = module_ci_metab_norm_args(),
                                        internal_standard_regex = "^IS_") {
  createMetaboliteAssayData(
    metabolite_data = module_ci_metab_norm_assays(
      layout = layout,
      samples = samples,
      include_itsd = include_itsd,
      include_controls = include_controls,
      positive_only = positive_only
    ),
    design_matrix = module_ci_metab_norm_design(samples = samples),
    metabolite_id_column = "database_identifier",
    annotation_id_column = "Name",
    sample_id = "Run",
    group_id = "group",
    technical_replicate_id = "tech_rep_group",
    database_identifier_type = "database_identifier",
    internal_standard_regex = internal_standard_regex,
    args = args
  )
}

module_ci_metab_norm_sample_cols <- function(object, assay_name = names(object@metabolite_data)[[1]]) {
  intersect(colnames(object@metabolite_data[[assay_name]]), as.character(object@design_matrix[[object@sample_id]]))
}

module_ci_metab_norm_matrix <- function(object, assay_name = names(object@metabolite_data)[[1]]) {
  assay <- object@metabolite_data[[assay_name]]
  sample_cols <- module_ci_metab_norm_sample_cols(object, assay_name)
  mat <- as.matrix(assay[, sample_cols, drop = FALSE])
  rownames(mat) <- as.character(assay[[object@metabolite_id_column]])
  mat
}

module_ci_metab_norm_feature_ids <- function(object, assay_name = names(object@metabolite_data)[[1]]) {
  as.character(object@metabolite_data[[assay_name]][[object@metabolite_id_column]])
}

module_ci_metab_norm_assert_alignment <- function(object, expected_assays = names(object@metabolite_data)) {
  testthat::expect_identical(names(object@metabolite_data), expected_assays)
  for (assay_name in expected_assays) {
    testthat::expect_setequal(
      module_ci_metab_norm_sample_cols(object, assay_name),
      as.character(object@design_matrix[[object@sample_id]])
    )
  }
  invisible(TRUE)
}

module_ci_metab_norm_ctrl <- function(object) {
  stats::setNames(
    lapply(names(object@metabolite_data), function(assay_name) {
      ids <- module_ci_metab_norm_feature_ids(object, assay_name)
      controls <- grepl("^CTRL_", ids)
      names(controls) <- ids
      controls
    }),
    names(object@metabolite_data)
  )
}

module_ci_metab_norm_counts <- function(object) {
  lapply(names(object@metabolite_data), function(assay_name) {
    assay <- object@metabolite_data[[assay_name]]
    list(
      features = nrow(assay),
      samples = length(module_ci_metab_norm_sample_cols(object, assay_name))
    )
  }) |> stats::setNames(names(object@metabolite_data))
}

module_ci_metab_norm_state_manager <- function(current_state,
                                               current_state_name = "metab_qc_complete") {
  states <- list()
  saved <- list()
  states[[current_state_name]] <- current_state

  manager <- new.env(parent = emptyenv())
  manager$current_state <- current_state_name
  manager$getState <- function(state = manager$current_state) {
    states[[state]]
  }
  manager$saveState <- function(state_name, s4_data_object, config_object, description) {
    states[[state_name]] <<- s4_data_object
    saved[[state_name]] <<- list(
      state_name = state_name,
      s4_data_object = s4_data_object,
      config_object = config_object,
      description = description
    )
    manager$current_state <- state_name
    invisible(TRUE)
  }
  manager$saved <- function() saved
  manager$states <- function() states
  manager
}

module_ci_metab_norm_workflow <- function(current_object = module_ci_metab_norm_object(),
                                          current_state_name = "metab_qc_complete") {
  list2env(list(
    state_manager = module_ci_metab_norm_state_manager(current_object, current_state_name),
    config_list = list(normalization = list(mode = "module-ci")),
    contrasts_tbl = data.frame(
      friendly_names = "KO_vs_WT",
      contrasts = "groupKO-groupWT",
      stringsAsFactors = FALSE
    ),
    design_matrix = current_object@design_matrix,
    metabolite_counts = module_ci_metab_norm_counts(current_object),
    qc_params = list(
      min_pearson_correlation_threshold = 0.75,
      internal_standard_regex = current_object@internal_standard_regex
    ),
    tab_status = list(
      setup_import = "complete",
      design_matrix = "complete",
      quality_control = "complete",
      normalization = "pending",
      differential_analysis = "locked"
    )
  ), parent = emptyenv())
}

module_ci_metab_norm_ruv_results <- function(object, mode = "manual") {
  ctrl <- module_ci_metab_norm_ctrl(object)
  stats::setNames(
    lapply(seq_along(ctrl), function(i) {
      assay_name <- names(ctrl)[[i]]
      list(
        success = TRUE,
        best_k = if (identical(mode, "automatic")) i else 2L,
        best_percentage = if (identical(mode, "automatic")) 12 + i else 40,
        best_realized_num_controls = sum(ctrl[[assay_name]], na.rm = TRUE),
        best_realized_percentage = 100 * sum(ctrl[[assay_name]], na.rm = TRUE) / nrow(object@metabolite_data[[assay_name]]),
        control_genes_index = ctrl[[assay_name]],
        cancor_plot = list(assay = assay_name, mode = mode),
        separation_score = if (identical(mode, "automatic")) 0.8 - (i / 100) else NA_real_,
        composite_score = if (identical(mode, "automatic")) 0.7 - (i / 100) else NA_real_,
        optimization_results = data.frame(
          percentage_requested = c(10, 11, 12),
          realized_num_controls = sum(ctrl[[assay_name]], na.rm = TRUE),
          best_k = c(3L, 2L, if (identical(mode, "automatic")) i else 2L),
          composite_score = c(0.3, 0.5, 0.7),
          status = "ok",
          stringsAsFactors = FALSE
        ),
        error = NULL
      )
    }),
    names(ctrl)
  )
}

module_ci_metab_norm_data <- function(object = module_ci_metab_norm_object()) {
  list2env(list(
    itsd_selections = lapply(object@metabolite_data, function(assay) {
      which(grepl("^IS_", assay[[object@metabolite_id_column]]))
    }),
    ruv_optimization_results = list(),
    correlation_results = list(),
    assay_names = names(object@metabolite_data),
    normalization_complete = FALSE,
    ruv_complete = FALSE,
    correlation_filtering_complete = FALSE,
    post_filter_obj = object,
    post_itsd_obj = NULL,
    post_log2_obj = NULL,
    post_norm_obj = NULL,
    ruv_corrected_obj = NULL,
    correlation_filtered_obj = NULL,
    plot_refresh_trigger = 0
  ), parent = emptyenv())
}

module_ci_metab_norm_input <- function(norm_method = "none",
                                       ruv_mode = "skip",
                                       apply_itsd = FALSE,
                                       itsd_aggregation = "median",
                                       log_offset = 1,
                                       correlation_threshold = 0.75,
                                       ruv_grouping_variable = "tech_rep_group") {
  list(
    norm_method = norm_method,
    ruv_mode = ruv_mode,
    apply_itsd = apply_itsd,
    itsd_aggregation = itsd_aggregation,
    log_offset = log_offset,
    min_pearson_correlation_threshold = correlation_threshold,
    ruv_grouping_variable = ruv_grouping_variable
  )
}

module_ci_metab_norm_correlation_pairs <- function(case = c("pass_all", "boundary", "fail_one", "fail_many", "small_n"),
                                                   assays = c("LCMS_Pos", "GCMS")) {
  case <- match.arg(case)
  correlations <- switch(
    case,
    pass_all = c(0.95, 0.90),
    boundary = c(0.75, 0.75),
    fail_one = c(0.95, 0.50),
    fail_many = c(0.40, 0.50),
    small_n = numeric(0)
  )

  stats::setNames(
    lapply(assays, function(assay_name) {
      if (length(correlations) == 0L) {
        return(data.frame(
          Run.x = character(),
          Run.y = character(),
          pearson_correlation = numeric(),
          stringsAsFactors = FALSE
        ))
      }
      data.frame(
        Run.x = c("S1", "S3"),
        Run.y = c("S2", "S4"),
        pearson_correlation = correlations,
        stringsAsFactors = FALSE
      )
    }),
    assays
  )
}

module_ci_metab_norm_export_variant <- function(variant = c("skip", "manual", "automatic")) {
  variant <- match.arg(variant)
  object <- suppressWarnings(normaliseUntransformedData(
    module_ci_metab_norm_object(positive_only = TRUE),
    method = "ITSD",
    itsd_aggregation = "median",
    remove_itsd_after_norm = TRUE
  ))
  object@args$daAnalysisParameters <- list(formula_string = "~ 0 + group + batch")

  workflow <- module_ci_metab_norm_workflow(
    current_object = object,
    current_state_name = if (identical(variant, "skip")) "metab_norm_complete" else "metab_correlation_filtered"
  )
  norm_data <- module_ci_metab_norm_data(object)
  norm_data$normalization_complete <- TRUE
  norm_data$ruv_complete <- TRUE
  norm_data$correlation_filtering_complete <- TRUE
  norm_data$post_norm_obj <- object
  norm_data$ruv_corrected_obj <- object
  norm_data$correlation_filtered_obj <- object
  norm_data$correlation_results <- module_ci_metab_norm_correlation_pairs("pass_all", names(object@metabolite_data))
  norm_data$ruv_optimization_results <- if (identical(variant, "skip")) {
    list()
  } else {
    module_ci_metab_norm_ruv_results(object, mode = variant)
  }

  buildMetabNormExportSessionData(
    workflowData = workflow,
    normData = norm_data,
    inputValues = module_ci_metab_norm_input(
      norm_method = "quantile",
      ruv_mode = variant,
      apply_itsd = TRUE,
      itsd_aggregation = "median"
    ),
    experimentLabel = paste("MCI-015", variant),
    exportTimestamp = as.POSIXct("2026-05-05 12:00:00", tz = "UTC")
  )
}
