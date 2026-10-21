test_that("audit column encoding preserves the exact historical representation", {
  reference <- function(column) {
    encoded <- ifelse(is.na(column), "<NA>", as.character(column))
    paste0(nchar(encoded, type = "bytes"), ":", encoded)
  }
  latin <- iconv("caf\u00e9", to = "latin1")
  bytes <- rawToChar(as.raw(c(99, 97, 102, 255)))
  Encoding(bytes) <- "bytes"
  cases <- list(character(), numeric(), integer(), logical(),
    rep(NA_character_, 3L), rep(NA_real_, 3L), c(bytes, NA_character_),
    c(NA_character_, "<NA>", "", "a\u001fb", latin, "caf\u00e9"),
    c(NA_real_, NaN, Inf, -Inf, -0, 1e-100, 1e100, pi),
    c(NA_integer_, -1L, 0L, 1L), c(TRUE, FALSE, NA), as.raw(c(0, 10, 255)),
    c(1 + 2i, NA_complex_), factor(c("b", NA, "a"), levels = c("b", "a", "unused")),
    as.Date(c("2026-01-01", NA)), as.POSIXct(c("2026-01-01", NA), tz = "UTC"),
    matrix(c(1, NA, 2, 3), 2), list(1:2, NULL, NA, "x"))
  for (column in cases) {
    before <- serialize(column, NULL)
    expect_identical(encodePeptideQcKeyColumn(column), reference(column))
    expect_identical(serialize(encodePeptideQcKeyColumn(column), NULL, version = 2),
      serialize(reference(column), NULL, version = 2))
    expect_identical(serialize(column, NULL), before)
  }
  expect_identical(encodePeptideQcKeyColumn(stats::setNames(c(NA, "x"), c("a", "b"))),
    reference(stats::setNames(c(NA, "x"), c("a", "b"))))
})

test_that("trial candidate uses one substantive digest for both immutable audit identifiers", {
  skip_if(identical(Sys.getenv("MULTISCHOLAR_DIA025_AUDIT_DIAGNOSTIC"), "true"),
    "unchanged-code diagnostic, not a candidate claim")
  skip_if(!identical(Sys.getenv("MULTISCHOLAR_DIA025_AUDIT_MEASURE"), "true"),
    "rejected audit candidate is retained only as a guarded experiment")
  before <- module_ci_prot_peptide_object()
  after <- before
  after@peptide_data <- module_ci_prot_qvalue_filter(before@peptide_data)
  digests <- 0L
  original <- .peptideQcDigest
  audited <- testthat::with_mocked_bindings({
    .emitPeptideQcAuditRecord(before, after, "qvalue_filter", "raw", "qvalue_filtered",
      now = as.POSIXct("2026-01-01", tz = "UTC"))
  }, .peptideQcDigest = function(value) {
    if (is.list(value) && identical(value$stage_id, "qvalue_filter")) digests <<- digests + 1L
    original(value)
  })
  expect_identical(digests, 1L)
  record <- audited@args$peptide_qc_audit$records[[1L]]
  substantive <- record[!names(record) %in% c("record_id", "timestamp_utc", "canonical_digest")]
  expect_identical(record$canonical_digest, original(substantive))
  expect_identical(record$record_id, paste0("pqc:", substr(record$canonical_digest, 1L, 24L)))
})

test_that("substantive digest reuse preserves the production Apply and reduces its allocation work", {
  receipt <- Sys.getenv("MULTISCHOLAR_DIA008_RECEIPT")
  diagnostic <- identical(Sys.getenv("MULTISCHOLAR_DIA025_AUDIT_DIAGNOSTIC"), "true")
  reference_library <- Sys.getenv("MULTISCHOLAR_DIA008_B0_LIBRARY")
  skip_if(!identical(Sys.getenv("MULTISCHOLAR_DIA025_AUDIT_MEASURE"), "true"),
    "isolated audit comparison is run once by its fixed guarded selector")
  skip_if(!nzchar(receipt) || !nzchar(reference_library), "guarded predecessor comparison only")
  frozen <- new.env(parent = baseenv())
  lazyLoad(file.path(reference_library, "MultiScholaR", "R", "MultiScholaR"), envir = frozen)
  reference <- get(".emitPeptideQcAuditRecord", frozen)
  reference_import <- get(".peptideQcImmutableImport", frozen)
  # Isolate this emitter change; all scientific and storage helpers stay matched.
  environment(reference) <- asNamespace("MultiScholaR")
  environment(reference_import) <- asNamespace("MultiScholaR")
  candidate <- .emitPeptideQcAuditRecord
  candidate_import <- .peptideQcImmutableImport
  object <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
  object@peptide_data <- object@peptide_data[rep(seq_len(nrow(object@peptide_data)), 12L), ]
  object@peptide_data$Precursor.Id <- paste0(object@peptide_data$Precursor.Id,
    "_audit_", seq_len(nrow(object@peptide_data)))
  object@args$srlQvalueProteotypicPeptideClean <- list(
    input_matrix_column_ids = c("Run", "Precursor.Id", "Protein.Group", "Protein.Ids",
      "Stripped.Sequence", "Modified.Sequence", "Precursor.Charge",
      "Precursor.Quantity", "Precursor.Normalised"),
    exclude_contaminants = FALSE, qvalue_threshold = 0.01,
    global_qvalue_threshold = 0.01, global_pg_qvalue_threshold = 0.01,
    choose_only_proteotypic_peptide = 1, confidence_provenance_mode = "diann_main_report")
  expect_equal(nrow(object@peptide_data), 55296L)
  now <- function() as.POSIXct("2026-01-01", tz = "UTC")
  apply <- function(manager) {
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- manager
    config <- new.env(parent = emptyenv())
    config$config_list <- object@args
    update <- function(...) updateConfigParameter(..., env = config)
    withDIAStateScope(manager, function() runPeptideQvalueApplyStep(workflow,
      0.01, 0.01, TRUE, updateConfigParameterFn = update, nowFn = now))
    # Exercise the downstream cold read, not only the inner digest calculation.
    manager$getState()
  }
  measure <- function(emitter, import_fingerprint, label) {
    manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(), "audit", "DIA"))
    on.exit(closeDIAStateAdapter(manager), add = TRUE)
    manager$saveState("raw_data_s4", object, list(), "source")
    log <- file.path(dirname(receipt), paste0("audit-", label, "-alloc.log"))
    calls <- c(emitter = 0L, substantive_digest = 0L, digest = 0L)
    actual_digest <- .peptideQcDigest
    initial <- gc(reset = TRUE)
    started <- proc.time()[["elapsed"]]
    utils::Rprofmem(log, threshold = 10000)
    on.exit(utils::Rprofmem(NULL), add = TRUE)
    result <- testthat::with_mocked_bindings(apply(manager),
      .emitPeptideQcAuditRecord = function(...) {
        calls[["emitter"]] <<- calls[["emitter"]] + 1L
        emitter(...)
      }, .peptideQcImmutableImport = import_fingerprint,
      .peptideQcDigest = function(value) {
        calls[["digest"]] <<- calls[["digest"]] + 1L
        if (is.list(value) && identical(value$stage_id, "qvalue_filter"))
          calls[["substantive_digest"]] <<- calls[["substantive_digest"]] + 1L
        actual_digest(value)
      })
    utils::Rprofmem(NULL)
    elapsed <- proc.time()[["elapsed"]] - started
    heap <- gc()
    expect_lt(file.size(log), 64 * 1024^2)
    allocation_lines <- grep("^[0-9]+ :", readLines(log), value = TRUE)
    sizes <- as.numeric(sub(" .*", "", allocation_lines))
    expect_identical(diaStateOwnership(manager)$working_records, 0L)
    expect_identical(manager$.context$store$open_resources, 0L)
    list(science_and_audit_sha256 = digest::digest(result, algo = "sha256"),
      calls = as.list(calls), audit_records = length(result@args$peptide_qc_audit$records),
      emitter_body_sha256 = digest::digest(deparse(body(emitter)), algo = "sha256"),
      largest_allocations = head(allocation_lines[order(sizes, decreasing = TRUE)], 12L),
      initial_vector_heap_mib = unname(initial["Vcells", 2L]),
      vector_heap_high_water_mib = unname(heap["Vcells", 6L]),
      allocations_over_10kb = sum(sizes), elapsed = elapsed,
      states = manager$getHistory())
  }
  # Warm the existing dispatch/storage path on the existing tiny fixture only.
  warm <- module_ci_prot_peptide_object()
  invisible(.emitPeptideQcAuditRecord(warm, warm, "qvalue_filter", "raw", "warm",
    status = "skipped", decision_reason = "warmup", now = now()))
  before <- measure(reference, reference_import, "before")
  after <- measure(candidate, candidate_import, "after")
  jsonlite::write_json(list(rows = nrow(object@peptide_data), samples = nrow(object@design_matrix),
    before = before, after = after, unchanged_code_diagnostic = diagnostic, performance_acceptance = FALSE,
    boundary = "production Apply helper, audit, durable checkpoint and downstream read; no GUI plots/startup"),
    file.path(dirname(receipt), "repair-audit-digest-reuse.json"),
    auto_unbox = TRUE, pretty = TRUE, digits = NA)
  expect_identical(after$science_and_audit_sha256, before$science_and_audit_sha256)
  expect_identical(after$states, before$states)
  expect_identical(before$calls$emitter, 1L)
  expect_identical(after$calls$emitter, 1L)
  expect_identical(before$audit_records, 1L)
  expect_identical(after$audit_records, 1L)
  if (diagnostic) {
    expect_identical(after$emitter_body_sha256, before$emitter_body_sha256)
    expect_identical(after$calls, before$calls)
  } else {
    expect_lt(after$allocations_over_10kb, before$allocations_over_10kb)
    expect_lte(after$vector_heap_high_water_mib, before$vector_heap_high_water_mib)
  }
})

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
  filtered <- before
  filtered@peptide_data <- module_ci_prot_qvalue_filter(before@peptide_data)
  manager$saveState("raw_data_s4", before, list(), "raw")
  save_filtered <- function(owner) .savePeptideQcState(
    owner, before, filtered, stage_id = "qvalue_filter",
    state_name = "qvalue_filtered", config_object = list(proteotypic_only = TRUE),
    description = "confidence filtered"
  )
  after <- testthat::with_mocked_bindings(save_filtered(manager),
    .emitPeptideQcAuditRecord = function(...) stop("audit must not execute"),
    .peptideQcDigest = function(...) stop("audit digest must not execute"))
  expect_identical(after, filtered)
  expect_null(after@args$peptide_qc_audit)
  expect_length(manager$getAuditRecords(), 0L)
  expect_identical(manager$getState(), filtered)
  expect_identical(manager$states$qvalue_filtered$config, list(proteotypic_only = TRUE))
  session_file <- withr::local_tempfile(fileext = ".rds")
  saveRDS(manager, session_file)
  restored <- readRDS(session_file)
  expect_false(restored$audit_enabled)
  expect_identical(restored$getState(), filtered)
  expect_identical(restored$revertToState("raw_data_s4"), before)

  opted_in <- WorkflowState$new(audit_enabled = TRUE)
  opted_in$saveState("raw_data_s4", before, list(), "raw")
  audited <- save_filtered(opted_in)
  expect_length(opted_in$getAuditRecords(), 1L)
  audited@args$peptide_qc_audit <- NULL
  expect_identical(audited, filtered)
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
  # KIT-046: the test-mode limpa object builder is gone; the protein object is built from the peptide lineage directly.
  quantity <- peptide@norm_quantity_column
  wide <- peptide_data |>
    dplyr::group_by(Protein.Group, Run) |>
    dplyr::summarise(value = sum(.data[[quantity]], na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = Run, values_from = value) |>
    as.data.frame(check.names = FALSE)
  protein <- ProteinQuantitativeData(
    protein_quant_table = wide,
    protein_id_column = "Protein.Group",
    protein_id_table = .proteinIdTableFromPeptideLineage(peptide, wide$Protein.Group),
    design_matrix = peptide@design_matrix,
    sample_id = peptide@sample_id,
    group_id = peptide@group_id,
    technical_replicate_id = peptide@technical_replicate_id,
    args = list()
  )
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
