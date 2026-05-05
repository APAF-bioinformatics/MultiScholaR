module_ci_prot_design_samples <- function(kind = "balanced") {
  switch(
    kind,
    sanitized = c("x123_sample_a", "wt_sample_2", "ko_sample_1", "ko_sample_2"),
    suffix = c("WT_01.raw", "WT_02.raw", "KO_01.raw", "KO_02.raw"),
    unbalanced = c("WT_1", "WT_2", "WT_3", "KO_1", "KO_2"),
    three_group = c("WT_1", "WT_2", "KO_1", "KO_2", "RES_1", "RES_2"),
    c("WT_1", "WT_2", "KO_1", "KO_2")
  )
}

module_ci_prot_design_data <- function(samples = module_ci_prot_design_samples()) {
  expand.grid(
    Protein.Ids = c("P001", "P002", "P003"),
    Run = samples,
    stringsAsFactors = FALSE
  ) |>
    transform(Intensity = seq_len(length(Protein.Ids)) * 100)
}

module_ci_prot_design_matrix <- function(kind = "two_group", samples = NULL) {
  if (is.null(samples)) {
    samples <- module_ci_prot_design_samples(if (identical(kind, "three_group")) "three_group" else kind)
  }
  group <- ifelse(grepl("^WT|wt|x123|WT_", samples), "WT", "KO")
  if (identical(kind, "three_group")) {
    group <- sub("_[0-9]+$", "", samples)
  }
  batch <- rep(c("B1", "B2"), length.out = length(samples))
  factor2 <- ifelse(batch == "B1", "low_batch", "high_batch")
  factor3 <- if (identical(kind, "interaction_ready")) {
    rep(c("female", "male"), length.out = length(samples))
  } else {
    NA_character_
  }

  design_group <- if (identical(kind, "multi_factor")) {
    paste(group, factor2, sep = "_")
  } else if (identical(kind, "interaction_ready")) {
    paste(group, factor2, factor3, sep = "_")
  } else {
    group
  }

  data.frame(
    Run = samples,
    factor1 = group,
    factor2 = factor2,
    factor3 = factor3,
    group = design_group,
    batch = batch,
    replicates = ave(seq_along(samples), group, FUN = seq_along),
    tech_reps = NA_integer_,
    stringsAsFactors = FALSE
  )
}

module_ci_prot_contrast_data <- function(...) {
  pairs <- list(...)
  if (length(pairs) == 0L) {
    pairs <- list(c("KO", "WT"))
  }
  do.call(rbind, lapply(pairs, function(pair) {
    data.frame(
      contrast_name = paste0(pair[[1]], ".vs.", pair[[2]]),
      numerator = pair[[1]],
      denominator = pair[[2]],
      stringsAsFactors = FALSE
    )
  }))
}

module_ci_prot_build_design_payload <- function(
    kind = "two_group",
    formula = "~ 0 + group",
    contrasts = module_ci_prot_contrast_data(c("KO", "WT")),
    removed = character()
) {
  design <- module_ci_prot_design_matrix(kind)
  data <- module_ci_prot_design_data(design$Run)
  config <- list(deAnalysisParameters = list(formula_string = "~ old"))
  buildProtDesignSaveResultsPayload(
    designMatrix = design,
    currentRemovedSamples = removed,
    dataCln = data,
    contrastData = contrasts,
    configList = config,
    formulaString = formula
  )
}

module_ci_prot_assert_design_samples_once <- function(design, expected_samples) {
  testthat::expect_setequal(as.character(design$Run), expected_samples)
  testthat::expect_false(anyDuplicated(as.character(design$Run)) > 0L)
  invisible(TRUE)
}

module_ci_prot_assert_design_rejected <- function(design, expected_samples) {
  actual <- as.character(design$Run)
  testthat::expect_true(
    anyDuplicated(actual) > 0L || !setequal(actual, expected_samples),
    info = "Invalid design should have duplicate, missing, or extra samples"
  )
  invisible(TRUE)
}

module_ci_prot_parse_contrast_terms <- function(contrast_string) {
  unlist(strsplit(gsub("\\s+", "", contrast_string), "[+-]"), use.names = FALSE)
}

module_ci_prot_assert_contrast_terms_in_model <- function(design, contrasts, formula = "~ 0 + group") {
  model_cols <- colnames(stats::model.matrix(stats::as.formula(formula), data = design))
  for (contrast in contrasts$contrasts) {
    terms <- module_ci_prot_parse_contrast_terms(contrast)
    testthat::expect_true(
      all(terms %in% model_cols),
      info = sprintf(
        "Contrast '%s' has term(s) outside model matrix columns [%s]",
        contrast,
        paste(model_cols, collapse = ", ")
      )
    )
  }
  invisible(model_cols)
}

module_ci_prot_validate_da_setup <- function(data, design, contrasts, formula = "~ 0 + group") {
  module_ci_prot_assert_design_samples_once(design, unique(as.character(data$Run)))
  module_ci_prot_assert_contrast_terms_in_model(design, contrasts, formula = formula)
  testthat::expect_true(all(c("Protein.Ids", "Run", "Intensity") %in% names(data)))
  testthat::expect_true(is.numeric(data$Intensity))
  invisible(TRUE)
}

module_ci_prot_design_workflow_data <- function(workflow_type = "DIA") {
  state_manager <- new.env(parent = emptyenv())
  state_manager$workflow_type <- workflow_type
  state_manager$saved <- list()
  state_manager$setWorkflowType <- function(type) {
    state_manager$workflow_type <- type
    invisible(type)
  }
  state_manager$saveState <- function(state_name, s4_data_object, config_object, description) {
    state_manager$saved[[state_name]] <- list(
      state_name = state_name,
      config_object = config_object,
      description = description
    )
    invisible(TRUE)
  }
  list2env(list(
    data_tbl = NULL,
    data_cln = NULL,
    design_matrix = NULL,
    contrasts_tbl = NULL,
    column_mapping = list(protein_col = "Protein.Ids", run_col = "Run", quantity_col = "Intensity"),
    config_list = list(
      globalParameters = list(workflow_type = workflow_type),
      deAnalysisParameters = list(formula_string = "~ 0 + group")
    ),
    taxon_id = 9606L,
    organism_name = "Homo sapiens",
    processing_log = list(),
    tab_status = list(
      setup_import = "complete",
      design_matrix = "pending",
      quality_control = "disabled",
      normalization = "disabled",
      differential_expression = "disabled",
      enrichment_analysis = "disabled",
      session_summary = "disabled"
    ),
    state_manager = state_manager
  ), parent = emptyenv())
}

module_ci_prot_write_design_import_pack <- function(root, payload) {
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  utils::write.table(payload$design_matrix, file.path(root, "design_matrix.tab"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(payload$data_cln, file.path(root, "data_cln.tab"), sep = "\t", quote = FALSE, row.names = FALSE)
  if (!is.null(payload$contrasts_tbl)) {
    writeLines(payload$contrasts_tbl$contrasts, file.path(root, "contrast_strings.tab"), useBytes = TRUE)
  }
  writeLines("[globalParameters]\nworkflow_type = DIA\n[deAnalysisParameters]\nformula_string = ~ 0 + group", file.path(root, "config.ini"))
  root
}

module_ci_prot_read_tsv <- function(path) {
  utils::read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
}
