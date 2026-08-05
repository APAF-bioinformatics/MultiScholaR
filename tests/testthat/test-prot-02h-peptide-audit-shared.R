test_that("peptide-QC audit records are deterministic apart from timestamps", {
  before_one <- module_ci_prot_peptide_object()
  before_two <- module_ci_prot_peptide_object()
  params <- list(
    qvalue_threshold = 0.01,
    global_qvalue_threshold = 0.01,
    global_pg_qvalue_threshold = 0.01,
    proteotypic_only = TRUE
  )
  filtered_one <- module_ci_prot_qvalue_filter(before_one@peptide_data)
  filtered_two <- module_ci_prot_qvalue_filter(before_two@peptide_data)
  after_one <- before_one
  after_two <- before_two
  after_one@peptide_data <- filtered_one
  after_two@peptide_data <- filtered_two

  audited_one <- .emitPeptideQcAuditRecord(
    before_one,
    after_one,
    stage_id = "qvalue_filter",
    parent_state = "raw_data_s4",
    current_state = "qvalue_filtered",
    resolved_parameters = params,
    now = as.POSIXct("2026-01-01 00:00:00", tz = "UTC")
  )
  audited_two <- .emitPeptideQcAuditRecord(
    before_two,
    after_two,
    stage_id = "qvalue_filter",
    parent_state = "raw_data_s4",
    current_state = "qvalue_filtered",
    resolved_parameters = params,
    now = as.POSIXct("2026-07-01 00:00:00", tz = "UTC")
  )

  record_one <- audited_one@args$peptide_qc_audit$records[[1L]]
  record_two <- audited_two@args$peptide_qc_audit$records[[1L]]
  expect_false(identical(record_one$timestamp_utc, record_two$timestamp_utc))
  expect_identical(record_one$record_id, record_two$record_id)
  expect_identical(record_one$canonical_digest, record_two$canonical_digest)
  expect_identical(
    .peptideQcSubstantiveAudit(audited_one),
    .peptideQcSubstantiveAudit(audited_two)
  )
  expect_identical(record_one$before_summary$rows, 20L)
  expect_identical(record_one$after_summary$rows, 8L)
  expect_identical(
    record_one$after_summary$frozen_identification_evidence$status,
    "frozen_after_identification_filter"
  )
  expect_equal(nrow(record_one$removal_ledger), 12L)
  expect_true(all(nzchar(record_one$removal_ledger$failure_reason)))
  expect_true(any(grepl("Q.Value_above_0.01", record_one$removal_ledger$failure_reason, fixed = TRUE)))
  expect_identical(
    record_one$immutable_import_digest,
    audited_one@args$peptide_qc_audit$immutable_import$digest
  )
})

test_that("confidence failure reasons report the configured thresholds", {
  data <- data.frame(
    Q.Value = c(0.051, 0.001),
    Global.Q.Value = c(0.001, 0.026),
    Global.PG.Q.Value = c(0.001, 0.001),
    Proteotypic = c(1, 1)
  )

  reasons <- .peptideQcConfidenceFailureReasons(
    data,
    list(
      qvalue_threshold = 0.05,
      global_qvalue_threshold = 0.025,
      global_pg_qvalue_threshold = 0.01,
      proteotypic_only = TRUE
    )
  )

  expect_identical(
    reasons,
    c("Q.Value_above_0.05", "Global.Q.Value_above_0.025")
  )
})

test_that("audit summaries keep frozen identification evidence distinct from survivors", {
  before <- module_ci_prot_peptide_object()
  annotated <- .annotateProteinIdentificationEvidence(
    before@peptide_data,
    protein_id_column = "Protein.Ids",
    peptide_sequence_column = "Stripped.Sequence",
    modified_peptide_sequence_column = "Modified.Sequence"
  )
  before@peptide_data <- annotated
  after <- before
  after@peptide_data <- after@peptide_data[
    after@peptide_data$Protein.Ids == "P_PASS",
    ,
    drop = FALSE
  ]

  audited <- .emitPeptideQcAuditRecord(
    before,
    after,
    stage_id = "protein_evidence_filter",
    parent_state = "precursor_rollup",
    current_state = "protein_peptide_filtered",
    resolved_parameters = list(
      num_peptides_per_protein_thresh = 1L,
      num_peptidoforms_per_protein_thresh = 2L
    )
  )
  summary <- tail(audited@args$peptide_qc_audit$records, 1L)[[1L]]$after_summary
  expect_identical(summary$frozen_identification_evidence$protein_count, 1L)
  expect_identical(summary$current_quantitative_survivors$active_protein_groups, 1L)
  expect_true("frozen_identification_evidence" %in% names(summary))
  expect_true("current_quantitative_survivors" %in% names(summary))
})

test_that("imputation audit ledger identifies every imputed cell", {
  before <- module_ci_prot_peptide_object()
  after <- before
  after@peptide_data$Peptide.Imputed <- after@peptide_data$Precursor.Normalised
  after@peptide_data$is_imputed <- FALSE
  target <- which(is.na(after@peptide_data$Peptide.Imputed))[[1L]]
  after@peptide_data$Peptide.Imputed[[target]] <- 42
  after@peptide_data$is_imputed[[target]] <- TRUE

  audited <- .emitPeptideQcAuditRecord(
    before,
    after,
    stage_id = "imputation",
    parent_state = "replicate_filtered",
    current_state = "imputed",
    resolved_parameters = list(proportion_missing_values = 0.75),
    transformation_type = "imputation"
  )
  record <- tail(audited@args$peptide_qc_audit$records, 1L)[[1L]]
  expect_identical(record$after_summary$imputed_cells, 1L)
  expect_equal(nrow(record$imputation_ledger), 1L)
  expect_identical(record$imputation_ledger$Peptide.Imputed, 42)
  expect_identical(record$imputation_ledger$source_stage, "imputation")
})

test_that("WorkflowState reverts active audit lineage without deleting immutable records", {
  manager <- WorkflowState$new()
  raw <- module_ci_prot_peptide_object()
  manager$saveState("raw_data_s4", raw, list(), "raw")

  first <- .emitPeptideQcAuditRecord(
    raw,
    raw,
    stage_id = "intensity_filter",
    parent_state = "raw_data_s4",
    current_state = "intensity_filtered",
    resolved_parameters = list(min_groups = 1L),
    status = "no_op"
  )
  manager$saveState("intensity_filtered", first, list(), "first")
  first_id <- first@args$peptide_qc_audit$current_record_id

  second <- .emitPeptideQcAuditRecord(
    first,
    first,
    stage_id = "sample_filter",
    parent_state = "intensity_filtered",
    current_state = "sample_filtered",
    resolved_parameters = list(min_peptides_per_sample = 1L),
    status = "no_op"
  )
  manager$saveState("sample_filtered", second, list(), "second")
  second_id <- second@args$peptide_qc_audit$current_record_id

  expect_setequal(names(manager$getAuditRecords()), c(first_id, second_id))
  manager$revertToState("intensity_filtered")
  expect_identical(manager$getStateAudit()$record_id, first_id)
  expect_setequal(names(manager$getAuditRecords()), c(first_id, second_id))
  expect_identical(manager$getStateAudit("raw_data_s4")$provenance_status, "legacy_or_not_applicable")
})

test_that("audit capture can be disabled without changing scientific data", {
  manager <- WorkflowState$new(audit_enabled = FALSE)
  before <- module_ci_prot_peptide_object()
  after <- .savePeptideQcState(
    manager,
    before,
    before,
    stage_id = "intensity_filter",
    state_name = "intensity_filtered",
    config_object = list(min_groups = 1L),
    description = "disabled audit"
  )
  expect_identical(after@peptide_data, before@peptide_data)
  expect_null(after@args$peptide_qc_audit)
  expect_length(manager$getAuditRecords(), 0L)
})

test_that("stage contracts fail with an actionable shape diagnostic", {
  object <- module_ci_prot_peptide_object()
  object@peptide_data$Global.PG.Q.Value <- NULL
  expect_error(
    .validatePeptideQcStagePrerequisites(object, "qvalue_filter"),
    "missing required column(s): Global.PG.Q.Value",
    fixed = TRUE
  )
})

test_that("protein-rollup audit retains Protein.Group as the active identity", {
  peptide_data <- module_ci_prot_peptide_table()
  peptide_data$Protein.Group <- paste0("GROUP_", peptide_data$Protein.Ids)
  peptide <- PeptideQuantitativeDataDiann(
    peptide_data = peptide_data,
    design_matrix = module_ci_prot_peptide_design(),
    sample_id = "Run",
    group_id = "group",
    technical_replicate_id = "replicate_group",
    args = list()
  )
  protein <- buildProtTestModeLimpaProteinObject(peptide)
  expect_identical(protein@protein_id_column, "Protein.Group")

  audited <- .emitPeptideQcAuditRecord(
    peptide,
    protein,
    stage_id = "protein_rollup",
    parent_state = "imputed",
    current_state = "protein_s4_created",
    resolved_parameters = list(rollup_method = "test_sum"),
    transformation_type = "aggregation"
  )
  record <- tail(audited@args$peptide_qc_audit$records, 1L)[[1L]]
  expect_identical(record$column_roles$active_protein_id, "Protein.Group")
  expect_identical(record$after_summary$active_protein_groups, 5L)
})
