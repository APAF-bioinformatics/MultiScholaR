module_ci_prot_peptide_design <- function(samples = c("A1", "A2", "B1", "B2")) {
  data.frame(
    Run = samples,
    group = sub("[0-9]+$", "", samples),
    replicate_group = sub("[0-9]+$", "", samples),
    batch = rep(c("B1", "B2"), length.out = length(samples)),
    stringsAsFactors = FALSE
  )
}

module_ci_prot_peptide_table <- function(values = NULL, q_values = NULL, global_q_values = NULL, proteotypic = NULL) {
  samples <- c("A1", "A2", "B1", "B2")
  peptides <- data.frame(
    Protein.Ids = c("P_PASS", "P_BOUNDARY", "P_MISSING", "P_ZERO", "P_INF"),
    Stripped.Sequence = c("PEP_PASS", "PEP_BOUNDARY", "PEP_MISSING", "PEP_ZERO", "PEP_INF"),
    stringsAsFactors = FALSE
  )

  if (is.null(values)) {
    values <- rbind(
      PEP_PASS = c(100, 110, 120, 130),
      PEP_BOUNDARY = c(10, 10, 10, 10),
      PEP_MISSING = c(NA, 20, NA, 25),
      PEP_ZERO = c(0, 0, 0, 0),
      PEP_INF = c(Inf, Inf, Inf, Inf)
    )
  }
  stopifnot(nrow(values) == nrow(peptides), ncol(values) == length(samples))

  rows <- expand.grid(
    peptide_index = seq_len(nrow(peptides)),
    Run = samples,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  rows$Protein.Ids <- peptides$Protein.Ids[rows$peptide_index]
  rows$Stripped.Sequence <- peptides$Stripped.Sequence[rows$peptide_index]
  rows$Precursor.Quantity <- as.numeric(values[cbind(rows$peptide_index, match(rows$Run, samples))])
  rows$Precursor.Normalised <- rows$Precursor.Quantity
  rows$Precursor.Id <- paste(rows$Protein.Ids, rows$Stripped.Sequence, rows$Run, sep = "_")
  rows$Modified.Sequence <- rows$Stripped.Sequence
  rows$Precursor.Charge <- 2L

  if (is.null(q_values)) {
    q_values <- c(0.001, 0.010, 0.020, NA_real_, 0.005)
  }
  if (is.null(global_q_values)) {
    global_q_values <- c(0.001, 0.010, 0.020, 0.005, NA_real_)
  }
  if (is.null(proteotypic)) {
    proteotypic <- c(1, 1, 1, 1, 0)
  }

  rows$Q.Value <- q_values[rows$peptide_index]
  rows$Global.Q.Value <- global_q_values[rows$peptide_index]
  rows$Proteotypic <- proteotypic[rows$peptide_index]
  rows[c(
    "Run",
    "Precursor.Id",
    "Protein.Ids",
    "Stripped.Sequence",
    "Modified.Sequence",
    "Precursor.Charge",
    "Precursor.Quantity",
    "Precursor.Normalised",
    "Q.Value",
    "Global.Q.Value",
    "Proteotypic"
  )]
}

module_ci_prot_peptide_object <- function(peptide_data = module_ci_prot_peptide_table(),
                                          design = module_ci_prot_peptide_design(unique(peptide_data$Run)),
                                          is_logged_data = TRUE) {
  obj <- methods::new(
    "PeptideQuantitativeData",
    peptide_data = peptide_data,
    protein_id_column = "Protein.Ids",
    peptide_sequence_column = "Stripped.Sequence",
    q_value_column = "Q.Value",
    global_q_value_column = "Global.Q.Value",
    proteotypic_peptide_sequence_column = "Proteotypic",
    raw_quantity_column = "Precursor.Quantity",
    norm_quantity_column = "Precursor.Normalised",
    is_logged_data = is_logged_data,
    design_matrix = design,
    sample_id = "Run",
    group_id = "group",
    technical_replicate_id = "replicate_group",
    args = list(
      peptideIntensityFiltering = list(
        groupwise_percentage_cutoff = 50,
        max_groups_percentage_cutoff = 50,
        peptides_intensity_cutoff_percentile = 1,
        core_utilisation = NA
      ),
      srlQvalueProteotypicPeptideClean = list(
        input_matrix_column_ids = c(
          "Run",
          "Precursor.Id",
          "Protein.Ids",
          "Stripped.Sequence",
          "Modified.Sequence",
          "Precursor.Charge",
          "Precursor.Quantity",
          "Precursor.Normalised"
        ),
        qvalue_threshold = 0.01,
        global_qvalue_threshold = 0.01,
        choose_only_proteotypic_peptide = 1
      )
    )
  )
  calcPeptideMatrix(obj)
}

module_ci_prot_peptide_ids <- function(tbl) {
  sort(unique(paste(tbl$Protein.Ids, tbl$Stripped.Sequence, sep = "%")))
}

module_ci_prot_assert_peptide_alignment <- function(tbl, design) {
  testthat::expect_setequal(unique(as.character(tbl$Run)), as.character(design$Run))
  testthat::expect_false(any(is.na(tbl$Protein.Ids) | tbl$Protein.Ids == ""))
  testthat::expect_false(any(is.na(tbl$Stripped.Sequence) | tbl$Stripped.Sequence == ""))
  testthat::expect_false(anyDuplicated(tbl[c("Protein.Ids", "Stripped.Sequence", "Run")]) > 0L)
  invisible(TRUE)
}

module_ci_prot_validate_peptide_design <- function(tbl, design) {
  peptide_samples <- unique(as.character(tbl$Run))
  design_samples <- as.character(design$Run)
  if (!setequal(peptide_samples, design_samples) || anyDuplicated(design_samples) > 0L) {
    stop("Peptide QC prerequisite error: peptide runs and design runs must match exactly once.", call. = FALSE)
  }
  invisible(TRUE)
}

module_ci_prot_intensity_filter <- function(tbl = module_ci_prot_peptide_table(),
                                            design = module_ci_prot_peptide_design(),
                                            threshold = 10,
                                            group_cutoff = 50,
                                            max_groups_cutoff = 50) {
  suppressMessages(peptideIntensityFilteringHelper(
    input_table = tbl,
    design_matrix = design,
    min_peptide_intensity_threshold = threshold,
    sample_id_column = "Run",
    grouping_variable = "group",
    groupwise_percentage_cutoff = group_cutoff,
    max_groups_percentage_cutoff = max_groups_cutoff,
    protein_id_column = Protein.Ids,
    peptide_sequence_column = Stripped.Sequence,
    peptide_quantity_column = Precursor.Quantity,
    core_utilisation = NA
  ))
}

module_ci_prot_qvalue_filter <- function(tbl = module_ci_prot_peptide_table(),
                                         q_threshold = 0.01,
                                         global_threshold = 0.01,
                                         proteotypic_only = 1) {
  suppressMessages(srlQvalueProteotypicPeptideCleanHelper(
    input_table = tbl,
    qvalue_threshold = q_threshold,
    global_qvalue_threshold = global_threshold,
    choose_only_proteotypic_peptide = proteotypic_only,
    input_matrix_column_ids = c(
      "Run",
      "Precursor.Id",
      "Protein.Ids",
      "Stripped.Sequence",
      "Modified.Sequence",
      "Precursor.Charge",
      "Precursor.Quantity",
      "Precursor.Normalised"
    ),
    protein_id_column = Protein.Ids,
    q_value_column = Q.Value,
    global_q_value_column = Global.Q.Value,
    proteotypic_peptide_sequence_column = Proteotypic
  ))
}

module_ci_prot_peptide_imputation_fixture <- function() {
  list(
    data = data.frame(
      Run = c("A1", "A2", "B1", "B2", "A1", "B1"),
      Protein.Ids = c("P1", "P1", "P1", "P1", "P2", "P3"),
      Stripped.Sequence = c("PEP1", "PEP1", "PEP1", "PEP1", "ALL_MISSING", "ONE_ONLY"),
      Peptide.Normalised = c(10, NA, 30, NA, NA, 5),
      stringsAsFactors = FALSE
    ),
    design = data.frame(
      Run = c("A1", "A2", "B1", "B2"),
      replicate_group = c("A", "A", "B", "B"),
      group = c("A", "A", "B", "B"),
      stringsAsFactors = FALSE
    )
  )
}

module_ci_prot_impute <- function(fixture = module_ci_prot_peptide_imputation_fixture(),
                                  proportion_missing_values = 0.75,
                                  core_utilisation = NA) {
  suppressMessages(peptideMissingValueImputationHelper(
    input_table = fixture$data,
    metadata_table = fixture$design,
    input_table_sample_id_column = Run,
    sample_id_tbl_sample_id_column = Run,
    replicate_group_column = replicate_group,
    protein_id_column = Protein.Ids,
    peptide_sequence_column = Stripped.Sequence,
    quantity_to_impute_column = Peptide.Normalised,
    imputed_value_column = Peptide.Imputed,
    hek_string = "HEK",
    proportion_missing_values = proportion_missing_values,
    core_utilisation = core_utilisation
  ))
}

module_ci_prot_state_manager <- function(initial_state) {
  state <- initial_state
  saved <- list()
  manager <- new.env(parent = emptyenv())
  manager$getState <- function(...) state
  manager$saveState <- function(state_name, s4_data_object, config_object, description) {
    saved[[length(saved) + 1L]] <<- list(
      state_name = state_name,
      s4_data_object = s4_data_object,
      config_object = config_object,
      description = description
    )
    state <<- s4_data_object
    invisible(TRUE)
  }
  manager$getHistory <- function() vapply(saved, `[[`, character(1), "state_name")
  manager$saved <- function() saved
  manager
}

module_ci_prot_workflow <- function(initial_state = module_ci_prot_peptide_object()) {
  workflow <- new.env(parent = emptyenv())
  workflow$state_manager <- module_ci_prot_state_manager(initial_state)
  workflow$qc_params <- list()
  workflow
}

module_ci_prot_update_config <- function(theObject, function_name, parameter_name, new_value) {
  if (is.null(theObject@args[[function_name]])) {
    theObject@args[[function_name]] <- list()
  }
  theObject@args[[function_name]][[parameter_name]] <- new_value
  theObject
}

module_ci_prot_update_missing <- function(theObject, min_reps_per_group, min_groups, function_name, grouping_variable) {
  if (is.null(theObject@args[[function_name]])) {
    theObject@args[[function_name]] <- list()
  }
  theObject@args[[function_name]]$groupwise_percentage_cutoff <- 50
  theObject@args[[function_name]]$max_groups_percentage_cutoff <- 50
  theObject
}

module_ci_prot_validate_peptide_prereqs <- function(workflow, required_state = "imputed") {
  history <- workflow$state_manager$getHistory()
  if (!required_state %in% history) {
    stop(
      sprintf("Protein QC prerequisite not met: peptide state '%s' is required before protein QC.", required_state),
      call. = FALSE
    )
  }
  invisible(TRUE)
}
