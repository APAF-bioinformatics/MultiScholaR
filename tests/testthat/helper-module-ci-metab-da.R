module_ci_metab_da_samples <- function() {
  paste0("S", seq_len(6))
}

module_ci_metab_da_design <- function(samples = module_ci_metab_da_samples()) {
  data.frame(
    Run = samples,
    group = rep(c("WT", "KO", "RES"), each = 2, length.out = length(samples)),
    batch = rep(c("B1", "B2"), length.out = length(samples)),
    replicates = rep(c("R1", "R2"), length.out = length(samples)),
    tech_rep_group = rep(c("Pair1", "Pair1", "Pair2", "Pair2", "Pair3", "Pair3"), length.out = length(samples)),
    stringsAsFactors = FALSE
  )
}

module_ci_metab_da_feature_ids <- function() {
  c("M_UP", "M_DOWN", "M_NS", "M_CONST", "M_TIED", "M_METALESS")
}

module_ci_metab_da_assay <- function(samples = module_ci_metab_da_samples(),
                                     assay_shift = 0,
                                     missing_metadata = FALSE) {
  values <- rbind(
    M_UP = c(10, 11, 31, 32, 18, 19),
    M_DOWN = c(50, 52, 22, 23, 43, 44),
    M_NS = c(20, 21, 22, 21, 20, 22),
    M_CONST = c(7, 7, 7, 7, 7, 7),
    M_TIED = c(13, 14, 28, 29, 23, 24),
    M_METALESS = c(16, 17, 12, 13, 19, 20)
  )
  values <- values[, seq_along(samples), drop = FALSE] + assay_shift
  colnames(values) <- samples

  sample_df <- as.data.frame(values, check.names = FALSE)
  names(sample_df) <- samples

  assay <- data.frame(
    database_identifier = rownames(values),
    stringsAsFactors = FALSE
  )
  if (!isTRUE(missing_metadata)) {
    assay$Name <- paste("Metabolite", rownames(values))
  }
  cbind(assay, sample_df)
}

module_ci_metab_da_assays <- function(layout = c("combined", "lc", "gc"),
                                      missing_metadata = FALSE) {
  layout <- match.arg(layout)
  lc <- module_ci_metab_da_assay(assay_shift = 0, missing_metadata = missing_metadata)
  gc <- module_ci_metab_da_assay(assay_shift = 5, missing_metadata = missing_metadata)

  switch(
    layout,
    lc = list(LCMS_Pos = lc),
    gc = list(GCMS = gc),
    combined = list(LCMS_Pos = lc, GCMS = gc)
  )
}

module_ci_metab_da_object <- function(layout = c("combined", "lc", "gc"),
                                      missing_metadata = FALSE,
                                      formula_string = "~ 0 + group") {
  layout <- match.arg(layout)
  methods::new(
    "MetaboliteAssayData",
    metabolite_data = module_ci_metab_da_assays(
      layout = layout,
      missing_metadata = missing_metadata
    ),
    metabolite_id_column = "database_identifier",
    annotation_id_column = if (isTRUE(missing_metadata)) "missing_annotation" else "Name",
    database_identifier_type = "database_identifier",
    internal_standard_regex = "^IS_",
    design_matrix = module_ci_metab_da_design(),
    sample_id = "Run",
    group_id = "group",
    technical_replicate_id = "tech_rep_group",
    args = list(
      daAnalysisParameters = list(formula_string = formula_string),
      module_ci_fixture = "MCI-016"
    )
  )
}

module_ci_metab_da_contrasts <- function(case = c(
  "two_group", "multi_group", "batch_aware", "invalid_term",
  "duplicate", "no_contrasts", "reversed", "contrast_string_column"
)) {
  case <- match.arg(case)
  switch(
    case,
    two_group = data.frame(
      friendly_names = "KO_vs_WT",
      contrasts = "groupKO-groupWT",
      full_format = "KO_vs_WT=groupKO-groupWT",
      stringsAsFactors = FALSE
    ),
    multi_group = data.frame(
      friendly_names = c("KO_vs_WT", "RES_vs_WT"),
      contrasts = c("groupKO-groupWT", "groupRES-groupWT"),
      full_format = c("KO_vs_WT=groupKO-groupWT", "RES_vs_WT=groupRES-groupWT"),
      stringsAsFactors = FALSE
    ),
    batch_aware = data.frame(
      friendly_names = "KO_vs_WT_batch_adjusted",
      contrasts = "groupKO-groupWT",
      full_format = "KO_vs_WT_batch_adjusted=groupKO-groupWT",
      stringsAsFactors = FALSE
    ),
    invalid_term = data.frame(
      friendly_names = "Missing_vs_WT",
      contrasts = "groupMISSING-groupWT",
      full_format = "Missing_vs_WT=groupMISSING-groupWT",
      stringsAsFactors = FALSE
    ),
    duplicate = data.frame(
      friendly_names = c("KO_vs_WT", "KO_vs_WT"),
      contrasts = c("groupKO-groupWT", "groupRES-groupWT"),
      full_format = c("KO_vs_WT=groupKO-groupWT", "KO_vs_WT=groupRES-groupWT"),
      stringsAsFactors = FALSE
    ),
    no_contrasts = data.frame(
      friendly_names = character(),
      contrasts = character(),
      full_format = character(),
      stringsAsFactors = FALSE
    ),
    reversed = data.frame(
      friendly_names = "WT_vs_KO",
      contrasts = "groupWT-groupKO",
      full_format = "WT_vs_KO=groupWT-groupKO",
      stringsAsFactors = FALSE
    ),
    contrast_string_column = data.frame(
      friendly_names = "KO_vs_WT",
      contrast_string = "groupKO-groupWT",
      stringsAsFactors = FALSE
    )
  )
}

module_ci_metab_da_state_manager <- function(current_state,
                                             current_state_name = "metab_correlation_filtered") {
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

module_ci_metab_da_workflow <- function(current_object = module_ci_metab_da_object(),
                                        contrasts = module_ci_metab_da_contrasts("multi_group")) {
  list2env(list(
    state_manager = module_ci_metab_da_state_manager(current_object),
    contrasts_tbl = contrasts,
    design_matrix = current_object@design_matrix,
    config_list = list(da = list(mode = "module-ci")),
    tab_status = list(
      setup_import = "complete",
      design_matrix = "complete",
      quality_control = "complete",
      normalization = "complete",
      differential_analysis = "pending"
    )
  ), parent = emptyenv())
}

module_ci_metab_da_data_env <- function() {
  list2env(list(
    current_s4_object = NULL,
    contrasts_available = NULL,
    assays_available = NULL,
    da_results_list = NULL,
    analysis_complete = FALSE,
    formula_from_s4 = NULL,
    current_row_clusters = NULL,
    current_col_clusters = NULL,
    current_heatmap_plot = NULL
  ), parent = emptyenv())
}

module_ci_metab_da_session_data <- function(layout = c("combined", "lc", "gc"),
                                            assay_names = NULL,
                                            missing_metadata = FALSE,
                                            contrasts = module_ci_metab_da_contrasts("multi_group")) {
  layout <- match.arg(layout)
  object <- module_ci_metab_da_object(layout = layout, missing_metadata = missing_metadata)
  list(
    omic_type = "metabolomics",
    current_s4_object = object,
    r6_current_state_name = "metab_correlation_filtered",
    contrasts_tbl = contrasts,
    assay_names = if (is.null(assay_names)) names(object@metabolite_data) else assay_names,
    feature_counts = vapply(object@metabolite_data, nrow, integer(1)),
    normalization_complete = TRUE,
    correlation_filtering_complete = TRUE
  )
}

module_ci_metab_da_capture_ui <- function() {
  capture <- new.env(parent = emptyenv())
  capture$select_updates <- list()
  capture$text_updates <- list()
  capture$notifications <- list()
  capture$removed <- character()
  capture$updateSelectInput <- function(session, inputId, choices, selected = NULL) {
    capture$select_updates[[length(capture$select_updates) + 1L]] <- list(
      session = session,
      inputId = inputId,
      choices = choices,
      selected = selected
    )
    invisible(NULL)
  }
  capture$updateTextAreaInput <- function(session, inputId, value) {
    capture$text_updates[[length(capture$text_updates) + 1L]] <- list(
      session = session,
      inputId = inputId,
      value = value
    )
    invisible(NULL)
  }
  capture$showNotification <- function(...) {
    capture$notifications[[length(capture$notifications) + 1L]] <- list(...)
    invisible(NULL)
  }
  capture$removeNotification <- function(id) {
    capture$removed <- c(capture$removed, id)
    invisible(NULL)
  }
  capture
}

module_ci_metab_da_local_binding <- function(env, name, value, .local_envir = parent.frame()) {
  had_binding <- exists(name, envir = env, inherits = FALSE)
  old_value <- if (had_binding) get(name, envir = env, inherits = FALSE) else NULL
  was_locked <- had_binding && bindingIsLocked(name, env)

  if (was_locked) {
    unlockBinding(name, env)
  }
  assign(name, value, envir = env)
  if (was_locked) {
    lockBinding(name, env)
  }

  withr::defer({
    if (exists(name, envir = env, inherits = FALSE) && bindingIsLocked(name, env)) {
      unlockBinding(name, env)
    }
    if (had_binding) {
      assign(name, old_value, envir = env)
    } else if (exists(name, envir = env, inherits = FALSE)) {
      rm(list = name, envir = env)
    }
    if (was_locked && exists(name, envir = env, inherits = FALSE)) {
      lockBinding(name, env)
    }
  }, envir = .local_envir)
}

module_ci_metab_da_stats_table <- function(feature_ids,
                                           case = c(
                                             "mixed", "no_significant", "all_significant",
                                             "no_variance", "tied_pvalues"
                                           )) {
  case <- match.arg(case)
  n <- length(feature_ids)
  log_fc <- stats::setNames(c(1.8, -1.6, 0.2, 0, 1.2, -0.9)[seq_len(n)], feature_ids)
  raw_p <- stats::setNames(c(0.001, 0.002, 0.40, 1, 0.006, 0.007)[seq_len(n)], feature_ids)
  fdr <- stats::setNames(c(0.01, 0.02, 0.50, 1, 0.03, 0.03)[seq_len(n)], feature_ids)

  if (identical(case, "no_significant")) {
    log_fc[] <- c(0.1, -0.1, 0.2, 0, 0.1, -0.1)[seq_len(n)]
    raw_p[] <- 0.60
    fdr[] <- 0.80
  } else if (identical(case, "all_significant")) {
    log_fc[] <- rep(c(1.5, -1.4), length.out = n)
    raw_p[] <- seq(0.001, by = 0.001, length.out = n)
    fdr[] <- seq(0.005, by = 0.001, length.out = n)
  } else if (identical(case, "no_variance")) {
    log_fc[] <- 0
    raw_p[] <- 1
    fdr[] <- 1
  } else if (identical(case, "tied_pvalues")) {
    log_fc[] <- rep(c(1.2, -1.2), length.out = n)
    raw_p[] <- 0.01
    fdr[] <- 0.02
  }

  data.frame(
    logFC = unname(log_fc),
    P.Value = unname(raw_p),
    raw_pvalue = unname(raw_p),
    fdr_value_bh = unname(fdr),
    fdr_qvalue = unname(fdr),
    row.names = feature_ids,
    check.names = FALSE
  )
}

module_ci_metab_da_bind_orchestration_stubs <- function(case = "mixed",
                                                        .local_envir = parent.frame()) {
  helper_env <- asNamespace("MultiScholaR")
  observed <- new.env(parent = emptyenv())
  observed$calls <- list()

  module_ci_metab_da_local_binding(
    helper_env,
    "runTestsContrastsMetabDA",
    function(data,
             contrast_strings,
             design_matrix,
             formula_string,
             sample_id_col = "Run",
             treat_lfc_cutoff = NA,
             eBayes_trend = TRUE,
             eBayes_robust = TRUE) {
      observed$calls[[length(observed$calls) + 1L]] <- list(
        sample_cols = colnames(data),
        contrast_strings = contrast_strings,
        formula_string = formula_string,
        sample_id_col = sample_id_col,
        treat_lfc_cutoff = treat_lfc_cutoff,
        eBayes_trend = eBayes_trend,
        eBayes_robust = eBayes_robust
      )
      stats <- stats::setNames(
        lapply(contrast_strings, function(...) {
          module_ci_metab_da_stats_table(rownames(data), case = case)
        }),
        contrast_strings
      )
      list(
        results = stats,
        fit.eb = structure(list(case = case), class = "module_ci_metab_da_fit"),
        qvalue_warnings = if (identical(case, "tied_pvalues")) {
          stats::setNames(as.list(rep(TRUE, length(contrast_strings))), contrast_strings)
        } else {
          list()
        }
      )
    },
    .local_envir = .local_envir
  )

  observed
}

module_ci_metab_da_results_long <- function(case = c(
  "mixed", "no_significant", "all_significant", "no_variance", "tied_pvalues"
)) {
  case <- match.arg(case)
  object <- module_ci_metab_da_object(layout = "combined")
  contrasts <- module_ci_metab_da_contrasts("multi_group")

  rows <- list()
  for (assay_name in names(object@metabolite_data)) {
    assay <- object@metabolite_data[[assay_name]]
    feature_ids <- as.character(assay$database_identifier)
    for (contrast_idx in seq_len(nrow(contrasts))) {
      stats <- module_ci_metab_da_stats_table(feature_ids, case = case)
      stats$metabolite_id <- rownames(stats)
      stats$metabolite_name <- paste("Metabolite", stats$metabolite_id)
      stats$assay <- assay_name
      stats$comparison <- contrasts$contrasts[[contrast_idx]]
      stats$friendly_name <- contrasts$friendly_names[[contrast_idx]]
      stats$significant <- ifelse(
        stats$fdr_qvalue < 0.05 & abs(stats$logFC) >= 0,
        ifelse(stats$logFC > 0, "Up", ifelse(stats$logFC < 0, "Down", "NS")),
        "NS"
      )
      stats$numerator <- sub("-.*$", "", stats$comparison)
      stats$denominator <- sub("^.*-", "", stats$comparison)
      for (sample in module_ci_metab_da_samples()) {
        group <- object@design_matrix$group[match(sample, object@design_matrix$Run)]
        col_name <- paste0("intensity.", sample, ".", group)
        stats[[col_name]] <- assay[match(stats$metabolite_id, assay$database_identifier), sample]
      }
      rows[[length(rows) + 1L]] <- stats
    }
  }

  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result[, c(
    "metabolite_id", "metabolite_name", "assay", "comparison",
    "friendly_name", "logFC", "raw_pvalue", "fdr_qvalue",
    "significant", "numerator", "denominator",
    grep("^intensity\\.", names(result), value = TRUE)
  )]
}

module_ci_metab_da_results_list <- function(case = "mixed") {
  long <- module_ci_metab_da_results_long(case = case)
  list(
    theObject = module_ci_metab_da_object(layout = "combined"),
    contrasts_results = list(),
    da_metabolites_long = long,
    per_assay_results = split(long, long$assay),
    significant_counts = lapply(split(long, long$assay), function(df) {
      list(
        up = sum(df$significant == "Up", na.rm = TRUE),
        down = sum(df$significant == "Down", na.rm = TRUE),
        ns = sum(df$significant == "NS", na.rm = TRUE)
      )
    }),
    qvalue_warnings = list()
  )
}

module_ci_metab_da_assert_long_schema <- function(results) {
  testthat::expect_s3_class(results, "data.frame")
  testthat::expect_true(all(c(
    "metabolite_id", "metabolite_name", "assay", "comparison",
    "friendly_name", "logFC", "raw_pvalue", "fdr_qvalue", "significant"
  ) %in% names(results)))
  testthat::expect_true(any(grepl("^intensity\\.", names(results))))
  testthat::expect_false(any(is.na(results$metabolite_id)))
  testthat::expect_false(any(is.na(results$assay)))
  testthat::expect_false(any(is.na(results$comparison)))
  invisible(TRUE)
}

module_ci_metab_da_report_schema_smoke <- function(paths) {
  tables <- lapply(paths, function(path) {
    table <- read.delim(path, check.names = FALSE)
    module_ci_metab_da_assert_long_schema(table)
    testthat::expect_true(all(nzchar(table$friendly_name)))
    testthat::expect_true(all(nzchar(table$comparison)))
    testthat::expect_true(all(grepl("^intensity\\.", names(table)) %in% c(TRUE, FALSE)))
    testthat::expect_true(any(grepl("^intensity\\.", names(table))))
    table
  })
  combined <- do.call(rbind, tables)
  list(
    table_count = length(tables),
    assays = sort(unique(combined$assay)),
    contrasts = sort(unique(combined$comparison)),
    intensity_columns = grep("^intensity\\.", names(combined), value = TRUE)
  )
}
