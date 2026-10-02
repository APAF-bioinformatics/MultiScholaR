# One member per guarded process: no sequential A/B heap-history confounding.
test_that("installed audit memory member completes science publication and readback", {
  receipt <- Sys.getenv("MULTISCHOLAR_DIA008_RECEIPT")
  skip_if(!identical(Sys.getenv("MULTISCHOLAR_DIA025_AUDIT_MEMBER"), "true"),
    "fixed guarded member selector only")
  library <- Sys.getenv("MULTISCHOLAR_DIA008_LIBRARY")
  expect_identical(normalizePath(getNamespaceInfo(asNamespace("MultiScholaR"), "path")),
    normalizePath(file.path(library, "MultiScholaR")))
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
  expect_identical(nrow(object@peptide_data), 55296L)
  input_sha256 <- digest::digest(object, algo = "sha256")
  manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(), "audit", "DIA"))
  on.exit(closeDIAStateAdapter(manager), add = TRUE)
  manager$saveState("raw_data_s4", object, list(), "source")
  workflow <- new.env(parent = emptyenv())
  workflow$state_manager <- manager
  config <- new.env(parent = emptyenv())
  config$config_list <- object@args
  update <- function(...) updateConfigParameter(..., env = config)
  now <- function() as.POSIXct("2026-01-01", tz = "UTC")
  started <- proc.time()[["elapsed"]]
  withDIAStateScope(manager, function() runPeptideQvalueApplyStep(workflow,
    0.01, 0.01, TRUE, updateConfigParameterFn = update, nowFn = now))
  result <- manager$getState()
  elapsed <- proc.time()[["elapsed"]] - started
  # Cumulative heap high-water includes startup/preparation. No pre-Apply GC or profiler.
  heap <- gc()
  status <- readLines("/proc/self/status")
  hwm <- as.numeric(sub("^VmHWM:[[:space:]]*([0-9]+).*", "\\1",
    grep("^VmHWM:", status, value = TRUE))) * 1024
  expect_length(result@args$peptide_qc_audit$records, 1L)
  expect_identical(diaStateOwnership(manager)$working_records, 0L)
  expect_identical(manager$.context$store$open_resources, 0L)
  report <- list(rows = nrow(object@peptide_data), samples = nrow(object@design_matrix),
    input_sha256 = input_sha256, science_and_audit_sha256 = digest::digest(result, algo = "sha256"),
    states = manager$getHistory(), apply_elapsed = elapsed,
    vector_heap_high_water_mib = unname(heap["Vcells", 6L]), main_r_hwm_bytes = hwm,
    profiling = FALSE, prewarm = FALSE, audit_records = 1L,
    boundary = "fresh installed process; production Apply, audit, checkpoint, readback; service includes startup/preparation",
    full_gui_acceptance = FALSE)
  jsonlite::write_json(report, file.path(dirname(receipt), "audit-member.json"),
    auto_unbox = TRUE, pretty = TRUE, digits = NA)
})
