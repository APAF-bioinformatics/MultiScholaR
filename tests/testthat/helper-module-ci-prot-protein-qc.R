module_ci_prot_protein_design <- function(samples = c("A1", "A2", "B1", "B2"),
                                          groups = c("A", "A", "B", "B"),
                                          replicates = c("R1", "R2", "R1", "R2")) {
  data.frame(
    Run = samples,
    group = groups,
    replicates = replicates,
    batch = rep(c("Batch1", "Batch2"), length.out = length(samples)),
    stringsAsFactors = FALSE
  )
}

module_ci_prot_protein_table <- function(proteins = c("P_PASS", "P_MISSING", "P_ZERO", "P_DUP", "P_DUP"),
                                         samples = c("A1", "A2", "B1", "B2"),
                                         values = NULL) {
  if (is.null(values)) {
    values <- rbind(
      P_PASS = c(20, 21, 22, 23),
      P_MISSING = c(NA, 25, NA, 26),
      P_ZERO = c(0, 0, 0, 0),
      P_DUP_1 = c(10, 12, 14, 16),
      P_DUP_2 = c(30, 32, 34, 36)
    )
  }
  stopifnot(nrow(values) == length(proteins), ncol(values) == length(samples))

  data.frame(
    Protein.Ids = proteins,
    as.data.frame(values, check.names = FALSE),
    check.names = FALSE,
    stringsAsFactors = FALSE
  ) |>
    stats::setNames(c("Protein.Ids", samples))
}

module_ci_prot_protein_object <- function(protein_quant_table = module_ci_prot_protein_table(),
                                          design = module_ci_prot_protein_design(colnames(protein_quant_table)[-1]),
                                          args = list()) {
  default_args <- list(
    removeRowsWithMissingValuesPercent = list(
      groupwise_percentage_cutoff = 50,
      max_groups_percentage_cutoff = 50,
      proteins_intensity_cutoff_percentile = 1
    ),
    ruvIII_C_Varying = list(
      ruv_grouping_variable = "group"
    ),
    limpa_dpc_quant_results = list(
      peptide_counts_per_protein = data.frame(
        Protein.Ids = protein_quant_table$Protein.Ids,
        peptide_count = rep(2L, nrow(protein_quant_table)),
        peptidoform_count = rep(2L, nrow(protein_quant_table)),
        stringsAsFactors = FALSE
      ),
      quantified_elist = NULL
    )
  )
  args <- utils::modifyList(default_args, args)

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

module_ci_prot_peptide_rollup_object <- function() {
  peptide_data <- data.frame(
    Run = rep(c("A1", "A2", "B1", "B2"), times = 4),
    Protein.Ids = rep(c("P_SHARED", "P_SHARED", "P_SINGLE", "P_AMBIG"), each = 4),
    Stripped.Sequence = rep(c("PEP_A", "PEP_B", "PEP_SINGLE", "PEP_AMBIG"), each = 4),
    Peptide.Imputed = c(
      10, 11, 12, 13,
      20, 21, 22, 23,
      5, NA, NA, NA,
      30, 31, 32, 33
    ),
    Peptide.Normalised = c(
      10, 11, 12, 13,
      20, 21, 22, 23,
      5, NA, NA, NA,
      30, 31, 32, 33
    ),
    Q.Value = rep(0.001, 16),
    PG.Q.Value = rep(0.001, 16),
    stringsAsFactors = FALSE
  )

  methods::new(
    "PeptideQuantitativeData",
    peptide_data = peptide_data,
    protein_id_column = "Protein.Ids",
    peptide_sequence_column = "Stripped.Sequence",
    q_value_column = "Q.Value",
    global_q_value_column = "PG.Q.Value",
    raw_quantity_column = "Peptide.Imputed",
    norm_quantity_column = "Peptide.Imputed",
    is_logged_data = TRUE,
    design_matrix = module_ci_prot_protein_design(),
    sample_id = "Run",
    group_id = "group",
    technical_replicate_id = "replicates",
    args = list(globalParameters = list(workflow_type = "DIA"))
  )
}

module_ci_prot_protein_state_manager <- function(initial_state = NULL,
                                                 initial_name = "initial") {
  state <- initial_state
  saved <- list()
  manager <- new.env(parent = emptyenv())
  manager$current_state <- initial_name
  manager$getState <- function(state_name = NULL) {
    if (!is.null(state_name)) {
      manager$current_state <<- state_name
    }
    state
  }
  manager$saveState <- function(state_name, s4_data_object, config_object, description) {
    saved[[length(saved) + 1L]] <<- list(
      state_name = state_name,
      s4_data_object = s4_data_object,
      config_object = config_object,
      description = description
    )
    state <<- s4_data_object
    manager$current_state <<- state_name
    invisible(TRUE)
  }
  manager$getHistory <- function() vapply(saved, `[[`, character(1), "state_name")
  manager$saved <- function() saved
  manager$revertToState <- function(state_name) {
    matched <- Filter(function(item) identical(item$state_name, state_name), saved)
    if (length(matched) == 0L) {
      stop(sprintf("State '%s' not found.", state_name), call. = FALSE)
    }
    state <<- matched[[length(matched)]]$s4_data_object
    manager$current_state <<- state_name
    state
  }
  manager
}

module_ci_prot_protein_workflow <- function(initial_state = module_ci_prot_protein_object(),
                                           initial_name = "protein_s4_created") {
  workflow <- new.env(parent = emptyenv())
  workflow$state_manager <- module_ci_prot_protein_state_manager(initial_state, initial_name)
  workflow$qc_params <- list()
  workflow$config_list <- list(globalParameters = list(workflow_type = "DIA"))
  workflow
}

module_ci_prot_protein_level_workflow <- function(workflow_type = "TMT",
                                                  data = module_ci_prot_protein_table(),
                                                  design = module_ci_prot_protein_design(colnames(data)[-1]),
                                                  protein_col = "Protein.Ids") {
  workflow <- new.env(parent = emptyenv())
  workflow$data_cln <- data
  workflow$design_matrix <- design
  workflow$column_mapping <- list(protein_col = protein_col)
  workflow$config_list <- list(globalParameters = list(workflow_type = workflow_type))
  workflow$state_manager <- module_ci_prot_protein_state_manager()
  workflow
}

module_ci_prot_iq_process_stub <- function(input_filename,
                                           output_filename,
                                           sample_id,
                                           primary_id,
                                           secondary_id,
                                           intensity_col,
                                           filter_double_less,
                                           normalization) {
  output <- data.frame(
    Protein.Ids = c("P_SHARED", "P_SINGLE"),
    S_001 = c(15, 5),
    S_002 = c(16, 0),
    S_003 = c(17, 0),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  readr::write_tsv(output, output_filename)
  invisible(output_filename)
}

module_ci_prot_accession_oracle <- function() {
  list(
    input = data.frame(
      row_id = 1:4,
      Protein.Ids = c(
        "P_LOW;P_BEST",
        "CON__BAD;P_ALT",
        "P_ISO-2;P_ISO",
        "P_NOGENE"
      ),
      stringsAsFactors = FALSE
    ),
    metadata = data.frame(
      uniprot_acc = c("P_LOW", "P_BEST", "P_ALT", "P_ISO", "P_NOGENE"),
      gene_name = c("LOW", "BEST", "ALT", "ISO", NA_character_),
      cleaned_acc = c("P_LOW", "P_BEST", "P_ALT", "P_ISO", "P_NOGENE"),
      protein_evidence = c(4, 1, 2, 1, 3),
      status = c("unreviewed", "reviewed", "reviewed", "reviewed", "unknown"),
      is_isoform = c("Canonical", "Canonical", "Canonical", "Canonical", "Canonical"),
      isoform_num = c(0, 0, 0, 0, 0),
      seq_length = c(100, 250, 150, 300, 50),
      annotation_score = c(1, 10, 5, 9, NA),
      stringsAsFactors = FALSE
    )
  )
}

module_ci_prot_assert_protein_alignment <- function(object) {
  sample_cols <- setdiff(colnames(object@protein_quant_table), object@protein_id_column)
  testthat::expect_setequal(sample_cols, as.character(object@design_matrix[[object@sample_id]]))
  testthat::expect_false(any(is.na(object@protein_quant_table[[object@protein_id_column]]) |
    object@protein_quant_table[[object@protein_id_column]] == ""))
  invisible(TRUE)
}

