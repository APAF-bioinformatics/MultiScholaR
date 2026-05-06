module_ci_lipid_qc_samples <- function() {
  c("S1", "S2", "S3", "S4")
}

module_ci_lipid_qc_design <- function(samples = module_ci_lipid_qc_samples(),
                                      missing_group = FALSE,
                                      small_n = FALSE) {
  groups <- if (isTRUE(small_n)) {
    c("WT", "KO", "RES", "REC")[seq_along(samples)]
  } else {
    rep(c("WT", "KO"), each = 2, length.out = length(samples))
  }
  if (isTRUE(missing_group)) {
    groups[[1]] <- NA_character_
  }
  data.frame(
    Run = samples,
    group = groups,
    batch = rep(c("B1", "B2"), length.out = length(samples)),
    tech_rep_group = paste0(groups, "_", rep(seq_len(2), each = 2, length.out = length(samples))),
    stringsAsFactors = FALSE
  )
}

module_ci_lipid_qc_assay <- function(samples = module_ci_lipid_qc_samples(),
                                     include_duplicates = TRUE,
                                     include_itsd = TRUE,
                                     constant = FALSE,
                                     nonfinite = FALSE,
                                     assay_shift = 0) {
  ids <- c("L_keep", "L_boundary", "L_zero", "L_missing")
  lipids <- c("Keep", "Boundary", "Zero", "Missing")
  classes <- c("PC", "TG", "SM", "PE")
  values <- rbind(
    c(100, 110, 120, 130),
    c(10, 10, 100, 100),
    c(0, 0, 0, 0),
    c(NA, 100, 100, 100)
  )

  if (isTRUE(include_duplicates)) {
    ids <- c(ids, "L_DUP", "L_DUP", "L_TIE", "L_TIE", "L_ANN_A", "L_ANN_B")
    lipids <- c(lipids, "Duplicate", "Duplicate", "Tie", "Tie", "SameAnnotation", "SameAnnotation")
    classes <- c(classes, "PC", "PC", "TG", "TG", "CE", "CE")
    values <- rbind(
      values,
      c(10, 15, 20, 25),
      c(90, 95, 100, 105),
      c(50, 50, 50, 50),
      c(50, 50, 50, 50),
      c(20, 22, 24, 26),
      c(80, 82, 84, 86)
    )
  }

  if (isTRUE(include_itsd)) {
    ids <- c(ids, "IS_PC_34_1", "IS_PC_34_1_ALT", "L_NAME_IS")
    lipids <- c(lipids, "IS_PC_34_1", "IS_PC_34_1", "IS_By_Name")
    classes <- c(classes, "ISTD", "ISTD", "ISTD")
    values <- rbind(
      values,
      c(1000, 1010, 990, 1005),
      c(900, 910, 890, 905),
      c(500, 520, 510, 530)
    )
  }

  if (isTRUE(nonfinite)) {
    ids <- c(ids, "L_INF", "L_NAN")
    lipids <- c(lipids, "Infinite", "NaN")
    classes <- c(classes, "PC", "TG")
    values <- rbind(
      values,
      c(Inf, 100, 100, 100),
      c(NaN, 100, 100, 100)
    )
  }

  if (isTRUE(constant)) {
    values[] <- 100
  }
  values <- values + assay_shift
  values[is.na(values)] <- NA
  sample_df <- as.data.frame(values, check.names = FALSE)
  names(sample_df) <- samples

  data.frame(
    lipid_id = ids,
    lipid = lipids,
    LipidClass = classes,
    sample_df,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

module_ci_lipid_qc_assays <- function(layout = "lc_pair",
                                      samples = module_ci_lipid_qc_samples(),
                                      include_duplicates = TRUE,
                                      include_itsd = TRUE,
                                      constant = FALSE) {
  primary <- module_ci_lipid_qc_assay(samples, include_duplicates, include_itsd, constant, assay_shift = 0)
  secondary <- module_ci_lipid_qc_assay(samples, include_duplicates, include_itsd, constant, assay_shift = 5)
  switch(
    layout,
    gc = list(GCMS = primary),
    combined = list(LCMS_Pos = primary, GCMS = secondary),
    lc_pair = list(LCMS_Pos = primary, LCMS_Neg = secondary),
    list(LCMS_Pos = primary)
  )
}

module_ci_lipid_qc_object <- function(layout = "lc_pair",
                                      samples = module_ci_lipid_qc_samples(),
                                      include_duplicates = TRUE,
                                      include_itsd = TRUE,
                                      constant = FALSE,
                                      missing_group = FALSE,
                                      small_n = FALSE,
                                      internal_standard_regex = "^IS_") {
  createLipidomicsAssayData(
    lipid_data = module_ci_lipid_qc_assays(
      layout = layout,
      samples = samples,
      include_duplicates = include_duplicates,
      include_itsd = include_itsd,
      constant = constant
    ),
    design_matrix = module_ci_lipid_qc_design(samples, missing_group = missing_group, small_n = small_n),
    lipid_id_column = "lipid_id",
    annotation_id_column = "lipid",
    sample_id = "Run",
    group_id = "group",
    technical_replicate_id = "tech_rep_group",
    database_identifier_type = "lipid_id",
    internal_standard_regex = internal_standard_regex,
    args = list()
  )
}

module_ci_lipid_qc_state_manager <- function(current_state) {
  saved <- list()
  history <- c("lipid_raw_data_s4", "lipid_intensity_filtered")
  manager <- new.env(parent = emptyenv())
  manager$getState <- function() current_state
  manager$saveState <- function(state_name, s4_data_object, config_object, description) {
    saved[[state_name]] <<- list(
      state_name = state_name,
      s4_data_object = s4_data_object,
      config_object = config_object,
      description = description
    )
    history <<- c(history, state_name)
    invisible(TRUE)
  }
  manager$getHistory <- function() history
  manager$saved <- function() saved
  manager
}

module_ci_lipid_qc_workflow <- function(current_state = module_ci_lipid_qc_object()) {
  list2env(list(
    state_manager = module_ci_lipid_qc_state_manager(current_state),
    config_list = list(qc = list(mode = "module-ci")),
    tab_status = list(
      setup_import = "complete",
      design_matrix = "complete",
      quality_control = "pending",
      normalization = "disabled"
    )
  ), parent = emptyenv())
}

module_ci_lipid_qc_feature_counts <- function(s4_object) {
  vapply(s4_object@lipid_data, nrow, integer(1))
}

module_ci_lipid_qc_ids <- function(assay) {
  as.character(assay$lipid_id)
}
