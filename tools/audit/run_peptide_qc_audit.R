#!/usr/bin/env Rscript

defaults <- list(
  input = NULL,
  output = NULL,
  design = NULL,
  config = NULL,
  contaminant_manifest = NULL,
  sample_column = "Run",
  group_column = "group",
  replicate_column = "replicates",
  full_pipeline = "false"
)

usage <- function() {
  cat(paste(
    "Usage: Rscript tools/audit/run_peptide_qc_audit.R --input <DIA-NN.tsv> --output <directory> [options]",
    "",
    "Options:",
    "  --design <TSV/CSV>                 Run/design metadata for design-dependent stages",
    "  --config <JSON/INI>                 Optional threshold overrides",
    "  --contaminant-manifest <TSV/CSV>    Optional contaminant accession manifest",
    "  --sample-column <name>              Input/design run column (default: Run)",
    "  --group-column <name>               Design biological group column (default: group)",
    "  --replicate-column <name>           Design technical replicate column (default: replicates)",
    "  --full-pipeline <true|false>         Apply design-dependent stages (requires --design)",
    "  --help",
    sep = "\n"
  ))
}

parse_args <- function(argv) {
  result <- defaults
  index <- 1L
  while (index <= length(argv)) {
    token <- argv[[index]]
    if (identical(token, "--help")) {
      usage()
      quit(status = 0L)
    }
    if (!startsWith(token, "--")) stop("Unexpected positional argument: ", token, call. = FALSE)
    split <- strsplit(sub("^--", "", token), "=", fixed = TRUE)[[1L]]
    key <- gsub("-", "_", split[[1L]], fixed = TRUE)
    if (!key %in% names(result)) stop("Unknown option: ", token, call. = FALSE)
    if (length(split) == 2L) {
      result[[key]] <- split[[2L]]
    } else {
      index <- index + 1L
      if (index > length(argv)) stop("Missing value for option: ", token, call. = FALSE)
      result[[key]] <- argv[[index]]
    }
    index <- index + 1L
  }
  result
}

as_flag <- function(value, name) {
  normalised <- tolower(trimws(as.character(value)))
  if (!normalised %in% c("true", "false", "1", "0", "yes", "no")) {
    stop(name, " must be true or false.", call. = FALSE)
  }
  normalised %in% c("true", "1", "yes")
}

read_table <- function(path) {
  delimiter <- if (grepl("\\.csv$", path, ignore.case = TRUE)) "," else "\t"
  as.data.frame(vroom::vroom(
    path,
    delim = delimiter,
    progress = FALSE,
    show_col_types = FALSE,
    .name_repair = "minimal"
  ), stringsAsFactors = FALSE)
}

read_config <- function(path) {
  if (is.null(path)) return(list())
  if (!file.exists(path)) stop("Config file does not exist.", call. = FALSE)
  if (grepl("\\.json$", path, ignore.case = TRUE)) {
    return(jsonlite::read_json(path, simplifyVector = TRUE))
  }
  ini::read.ini(path)
}

find_config_value <- function(config, names, default) {
  if (!is.list(config)) return(default)
  direct <- intersect(names, names(config))
  if (length(direct) > 0L) return(config[[direct[[1L]]]])
  for (value in config) {
    if (is.list(value)) {
      found <- find_config_value(value, names, structure(default, class = "pqc_default"))
      if (!inherits(found, "pqc_default")) return(found)
    }
  }
  if (inherits(default, "pqc_default")) default else default
}

as_positive_integer <- function(value, name) {
  value <- suppressWarnings(as.numeric(value))
  if (length(value) != 1L || is.na(value) || !is.finite(value) || value < 1 || value != floor(value)) {
    stop(name, " must be a positive integer.", call. = FALSE)
  }
  as.integer(value)
}

as_percentage <- function(value, name, include_zero = FALSE) {
  value <- suppressWarnings(as.numeric(value))
  lower_ok <- if (include_zero) value >= 0 else value > 0
  if (length(value) != 1L || is.na(value) || !is.finite(value) || !lower_ok || value > 100) {
    stop(name, " must be within the valid percentage range.", call. = FALSE)
  }
  value
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (is.null(args$input) || is.null(args$output)) {
  usage()
  stop("Both --input and --output are required.", call. = FALSE)
}
if (!file.exists(args$input)) stop("Input file does not exist.", call. = FALSE)
if (!file.exists("DESCRIPTION")) {
  stop("Run this command from the MultiScholaR repository root.", call. = FALSE)
}
full_pipeline <- as_flag(args$full_pipeline, "--full-pipeline")
if (full_pipeline && is.null(args$design)) {
  stop("--full-pipeline true requires --design.", call. = FALSE)
}
if (!is.null(args$design) && !file.exists(args$design)) stop("Design file does not exist.", call. = FALSE)
if (!is.null(args$contaminant_manifest) && !file.exists(args$contaminant_manifest)) {
  stop("Contaminant manifest does not exist.", call. = FALSE)
}

suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
.peptideQcDigest <- getFromNamespace(".peptideQcDigest", "MultiScholaR")
.peptideQcAuditSchemaVersion <- getFromNamespace(".peptideQcAuditSchemaVersion", "MultiScholaR")
.savePeptideQcState <- getFromNamespace(".savePeptideQcState", "MultiScholaR")
.peptideQcSubstantiveAudit <- getFromNamespace(".peptideQcSubstantiveAudit", "MultiScholaR")
`%||%` <- function(x, y) if (is.null(x)) y else x

raw_data <- read_table(args$input)
required_input <- c(
  args$sample_column,
  "Stripped.Sequence",
  "Modified.Sequence",
  "Q.Value",
  "Global.Q.Value",
  "Global.PG.Q.Value",
  "Proteotypic",
  "Precursor.Quantity"
)
if (!any(c("Protein.Group", "Protein.Ids") %in% names(raw_data))) {
  stop("Input requires Protein.Group or Protein.Ids.", call. = FALSE)
}
missing_input <- setdiff(required_input, names(raw_data))
if (length(missing_input) > 0L) {
  stop("Input is missing required DIA-NN column(s): ", paste(missing_input, collapse = ", "), call. = FALSE)
}
if (!"Precursor.Normalised" %in% names(raw_data)) {
  raw_data$Precursor.Normalised <- raw_data$Precursor.Quantity
}

design <- if (!is.null(args$design)) {
  read_table(args$design)
} else {
  runs <- unique(as.character(raw_data[[args$sample_column]]))
  data.frame(
    run = runs,
    group = "design_not_supplied",
    replicate = "design_not_supplied",
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) |>
    stats::setNames(c(args$sample_column, args$group_column, args$replicate_column))
}
required_design <- c(args$sample_column, args$group_column, args$replicate_column)
missing_design <- setdiff(required_design, names(design))
if (length(missing_design) > 0L) {
  stop("Design is missing required column(s): ", paste(missing_design, collapse = ", "), call. = FALSE)
}
if (anyDuplicated(as.character(design[[args$sample_column]])) > 0L) {
  stop("Design must contain exactly one row per run.", call. = FALSE)
}
if (!setequal(as.character(design[[args$sample_column]]), as.character(unique(raw_data[[args$sample_column]])))) {
  stop("Input and design run identities do not match exactly.", call. = FALSE)
}

config <- read_config(args$config)
parameters <- list(
  min_reps_per_group = as_positive_integer(
    find_config_value(config, c("min_reps_per_group"), 1L),
    "min_reps_per_group"
  ),
  min_groups = as_positive_integer(
    find_config_value(config, c("min_groups"), 1L),
    "min_groups"
  ),
  intensity_cutoff_percentile = as_percentage(
    find_config_value(config, c("peptides_intensity_cutoff_percentile", "intensity_cutoff_percentile"), 1),
    "intensity_cutoff_percentile",
    include_zero = TRUE
  ),
  min_peptides_per_protein = as_positive_integer(
    find_config_value(config, c("peptides_per_protein_cutoff", "num_peptides_per_protein_thresh"), 1L),
    "peptides_per_protein"
  ),
  min_peptidoforms_per_protein = as_positive_integer(
    find_config_value(config, c("peptidoforms_per_protein_cutoff", "num_peptidoforms_per_protein_thresh"), 2L),
    "peptidoforms_per_protein"
  ),
  min_peptides_per_sample = as_positive_integer(
    find_config_value(config, c("peptides_per_sample_cutoff"), 1L),
    "peptides_per_sample"
  ),
  imputation_missing_fraction = suppressWarnings(as.numeric(
    find_config_value(config, c("proportion_missing_values"), 0.75)
  ))
)
if (length(parameters$imputation_missing_fraction) != 1L ||
    is.na(parameters$imputation_missing_fraction) ||
    parameters$imputation_missing_fraction < 0 ||
    parameters$imputation_missing_fraction > 1) {
  stop("proportion_missing_values must be between 0 and 1.", call. = FALSE)
}

input_file_digest <- digest::digest(args$input, algo = "sha256", file = TRUE, serialize = FALSE)
parsed_import_digest <- .peptideQcDigest(raw_data)
active_protein_column <- if ("Protein.Group" %in% names(raw_data)) "Protein.Group" else "Protein.Ids"
q_args <- list(
  input_matrix_column_ids = names(raw_data),
  qvalue_threshold = 0.01,
  global_qvalue_threshold = 0.01,
  global_pg_qvalue_threshold = 0.01,
  choose_only_proteotypic_peptide = 1,
  confidence_provenance_mode = "diann_main_report",
  exclude_decoys = TRUE,
  exclude_contaminants = TRUE,
  contaminant_manifest = args$contaminant_manifest %||% ""
)
peptide <- methods::new(
  "PeptideQuantitativeData",
  peptide_data = raw_data,
  protein_id_column = active_protein_column,
  peptide_sequence_column = "Stripped.Sequence",
  q_value_column = "Q.Value",
  global_q_value_column = "Global.Q.Value",
  proteotypic_peptide_sequence_column = "Proteotypic",
  raw_quantity_column = "Precursor.Quantity",
  norm_quantity_column = "Precursor.Normalised",
  is_logged_data = FALSE,
  design_matrix = design,
  sample_id = args$sample_column,
  group_id = args$group_column,
  technical_replicate_id = args$replicate_column,
  args = list(srlQvalueProteotypicPeptideClean = q_args)
)
peptide@args$peptide_qc_audit <- list(
  schema_version = .peptideQcAuditSchemaVersion,
  provenance_status = "immutable_import",
  immutable_import = list(
    schema_version = .peptideQcAuditSchemaVersion,
    digest = input_file_digest,
    parsed_data_digest = parsed_import_digest,
    source_label = basename(args$input),
    source_path_recorded = FALSE
  ),
  records = list()
)

# The publication-facing identification baseline is deliberately computed
# before optional decoy/contaminant exclusion, then reported separately from
# the biological analysis view.
pre_biological <- suppressMessages(srlQvalueProteotypicPeptideClean(
  peptide,
  qvalue_threshold = 0.01,
  global_qvalue_threshold = 0.01,
  global_pg_qvalue_threshold = 0.01,
  choose_only_proteotypic_peptide = 1,
  input_matrix_column_ids = names(raw_data),
  confidence_provenance_mode = "diann_main_report",
  exclude_decoys = FALSE,
  exclude_contaminants = FALSE,
  contaminant_manifest = args$contaminant_manifest %||% ""
))
pre_biological_passing <- suppressMessages(filterMinNumPeptidesPerProtein(
  pre_biological,
  num_peptides_per_protein_thresh = parameters$min_peptides_per_protein,
  num_peptidoforms_per_protein_thresh = parameters$min_peptidoforms_per_protein
))

state_manager <- WorkflowState$new()
state_manager$saveState(
  "raw_data_s4",
  peptide,
  list(import_digest = input_file_digest, source_label = basename(args$input)),
  "Immutable DIA-NN import"
)
workflow <- new.env(parent = emptyenv())
workflow$state_manager <- state_manager
workflow$qc_params <- list()

no_log <- function(...) invisible(NULL)
update_config <- function(theObject, function_name, parameter_name, new_value) {
  if (!is.list(theObject@args[[function_name]])) theObject@args[[function_name]] <- list()
  theObject@args[[function_name]][[parameter_name]] <- new_value
  theObject
}
qvalue_result <- runPeptideQvalueApplyStep(
  workflow,
  qvalueThreshold = 0.01,
  globalQvalueThreshold = 0.01,
  globalPGQvalueThreshold = 0.01,
  proteotypicOnly = TRUE,
  updateConfigParameterFn = update_config,
  logInfoFn = no_log,
  logWarnFn = no_log
)
rollup_result <- runPeptideRollupApplyStep(workflow, logInfoFn = no_log)

save_skipped <- function(stage_id, state_name, reason, parameters = list()) {
  current <- state_manager$getState()
  saved <- .savePeptideQcState(
    state_manager = state_manager,
    before = current,
    after = current,
    stage_id = stage_id,
    state_name = state_name,
    config_object = parameters,
    description = paste("Skipped", stage_id, ":", reason),
    status = "skipped",
    decision_reason = reason,
    transformation_type = "no_op"
  )
  invisible(saved)
}

if (full_pipeline) {
  intensity_result <- runPeptideIntensityApplyStep(
    workflow,
    useStrictMode = FALSE,
    minRepsPerGroup = parameters$min_reps_per_group,
    minGroups = parameters$min_groups,
    intensityCutoffPercentile = parameters$intensity_cutoff_percentile,
    updateConfigParameterFn = update_config,
    logInfoFn = no_log
  )
} else {
  save_skipped(
    "intensity_filter",
    "intensity_filtered",
    "design-dependent quantitative filtering disabled; supply --design and --full-pipeline true",
    parameters[c("min_reps_per_group", "min_groups", "intensity_cutoff_percentile")]
  )
}

protein_evidence_result <- runProteinPeptideApplyStep(
  workflow,
  minPeptidesPerProtein = parameters$min_peptides_per_protein,
  minPeptidoformsPerProtein = parameters$min_peptidoforms_per_protein,
  updateConfigParameterFn = update_config,
  logInfoFn = no_log
)

if (full_pipeline) {
  sample_result <- runPeptideSampleApplyStep(
    workflow,
    minPeptidesPerSample = parameters$min_peptides_per_sample,
    updateConfigParameterFn = update_config,
    logInfoFn = no_log
  )
  replicate_result <- runPeptideReplicateApplyStep(
    workflow,
    replicateGroupColumn = args$replicate_column,
    logInfoFn = no_log
  )
  imputation_result <- runPeptideImputationStep(
    workflow,
    proportionMissingValues = parameters$imputation_missing_fraction,
    updateConfigParameterFn = update_config
  )
} else {
  save_skipped(
    "sample_filter", "sample_filtered",
    "design-dependent sample filtering disabled; supply --design and --full-pipeline true",
    list(min_peptides_per_sample = parameters$min_peptides_per_sample)
  )
  save_skipped(
    "replicate_filter", "replicate_filtered",
    "technical-replicate metadata was not supplied for execution",
    list(replicate_group_column = args$replicate_column, minimum_distinct_runs = 2L)
  )
  save_skipped(
    "imputation", "imputed",
    "technical-replicate imputation disabled without an executed design-dependent pipeline",
    list(proportion_missing_values = parameters$imputation_missing_fraction)
  )
}
save_skipped(
  "protein_rollup", "protein_s4_created",
  "headless identification audit does not execute an external protein rollup engine",
  list(active_protein_id = state_manager$getState()@protein_id_column)
)

final_object <- state_manager$getState()
audit <- .peptideQcSubstantiveAudit(final_object)
records <- audit$records
record_by_stage <- function(stage) {
  candidates <- Filter(function(record) identical(record$stage_id, stage), records)
  if (length(candidates) == 0L) NULL else tail(candidates, 1L)[[1L]]
}
q_record <- record_by_stage("qvalue_filter")
protein_record <- record_by_stage("protein_evidence_filter")
q_exclusion <- state_manager$getState("qvalue_filtered")@args$srlQvalueProteotypicPeptideClean$biological_exclusion_summary

pinned <- list(
  "b3bd682e1e34e7bb494be8d4162d25abb9cd6162f8f149818f6fd0e5ed0fba0a" = list(
    dataset = "cotton_report.tsv",
    three_gate_q_valid_proteotypic_rows = 18662L,
    three_gate_protein_groups = 1482L,
    passing_protein_groups_at_1_peptide_2_peptidoforms = 519L,
    legacy_two_gate_reference = list(rows = 19551L, protein_groups = 1673L)
  ),
  "d18f74c8a98e8ec9ec4ac0e9a4b32b6a1a0719e2b9ed46125e9e3f2aa267698e" = list(
    dataset = "KV_DIANN_report.tsv",
    three_gate_q_valid_proteotypic_rows = 159658L,
    three_gate_protein_groups = 2160L,
    passing_protein_groups_at_1_peptide_2_peptidoforms = 1763L,
    legacy_two_gate_reference = list(rows = 160121L, protein_groups = 2311L)
  )
)
observed_identification <- list(
  q_valid_proteotypic_rows_pre_biological_exclusion = nrow(pre_biological@peptide_data),
  protein_groups_pre_biological_exclusion = dplyr::n_distinct(pre_biological@peptide_data[[active_protein_column]]),
  passing_protein_groups_pre_biological_exclusion = dplyr::n_distinct(pre_biological_passing@peptide_data[[active_protein_column]]),
  q_valid_proteotypic_rows_after_biological_exclusion = q_record$after_summary$rows,
  protein_groups_after_biological_exclusion = q_record$after_summary$active_protein_groups,
  passing_protein_groups_after_biological_exclusion = protein_record$after_summary$active_protein_groups,
  biological_exclusion = if (is.data.frame(q_exclusion) && nrow(q_exclusion) > 0L) as.list(q_exclusion[1L, , drop = FALSE]) else list()
)
expected <- pinned[[input_file_digest]]
pinned_result <- if (is.null(expected)) {
  list(status = "input_hash_not_pinned")
} else {
  checks <- c(
    q_valid_proteotypic_rows = identical(
      as.integer(observed_identification$q_valid_proteotypic_rows_pre_biological_exclusion),
      as.integer(expected$three_gate_q_valid_proteotypic_rows)
    ),
    protein_groups = identical(
      as.integer(observed_identification$protein_groups_pre_biological_exclusion),
      as.integer(expected$three_gate_protein_groups)
    ),
    passing_protein_groups = identical(
      as.integer(observed_identification$passing_protein_groups_pre_biological_exclusion),
      as.integer(expected$passing_protein_groups_at_1_peptide_2_peptidoforms)
    )
  )
  list(status = if (all(checks)) "matched" else "mismatch", expected = expected, checks = as.list(checks))
}

summary <- list(
  schema_version = .peptideQcAuditSchemaVersion,
  source_label = basename(args$input),
  source_path_recorded = FALSE,
  input_sha256 = input_file_digest,
  parsed_import_digest = parsed_import_digest,
  immutable_import_digest_after_pipeline = audit$immutable_import$digest,
  immutable_import_verified = identical(input_file_digest, audit$immutable_import$digest) &&
    identical(parsed_import_digest, .peptideQcDigest(raw_data)),
  design_source = if (is.null(args$design)) "not_supplied" else "user_supplied",
  full_pipeline = full_pipeline,
  active_protein_id = active_protein_column,
  resolved_parameters = parameters,
  identification_closeout = observed_identification,
  pinned_baseline = pinned_result,
  state_history = state_manager$getHistory(),
  audit_record_ids = vapply(records, `[[`, character(1), "record_id"),
  canonical_audit_digest = .peptideQcDigest(audit)
)

dir.create(args$output, recursive = TRUE, showWarnings = FALSE)
saveRDS(audit, file.path(args$output, "peptide-qc-audit.rds"), version = 2)
jsonlite::write_json(
  summary,
  file.path(args$output, "peptide-qc-summary.json"),
  auto_unbox = TRUE,
  pretty = TRUE,
  na = "null",
  null = "null",
  digits = NA
)
cat(jsonlite::toJSON(summary, auto_unbox = TRUE, pretty = TRUE, na = "null", null = "null", digits = NA), "\n")
if (identical(pinned_result$status, "mismatch")) quit(status = 2L)
