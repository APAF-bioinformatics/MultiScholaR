# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

.peptideQcAuditSchemaVersion <- "1.0.0"

.peptideQcOr <- function(value, default) {
  if (is.null(value) || length(value) == 0L) default else value
}

.peptideQcCanonicalise <- function(value) {
  if (is.environment(value) || is.function(value) || identical(typeof(value), "externalptr")) {
    return(sprintf("<%s>", typeof(value)))
  }
  if (inherits(value, "POSIXt")) {
    return(format(value, "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC"))
  }
  if (is.factor(value)) {
    return(as.character(value))
  }
  if (is.data.frame(value)) {
    value <- as.data.frame(value, stringsAsFactors = FALSE)
    value <- value[, sort(names(value)), drop = FALSE]
    value[] <- lapply(value, .peptideQcCanonicalise)
    if (nrow(value) > 1L && ncol(value) > 0L) {
      row_keys <- do.call(
        paste,
        c(lapply(value, function(column) {
          encoded <- ifelse(is.na(column), "<NA>", as.character(column))
          paste0(nchar(encoded, type = "bytes"), ":", encoded)
        }), sep = "\u001f")
      )
      value <- value[order(row_keys, method = "radix"), , drop = FALSE]
    }
    rownames(value) <- NULL
    return(value)
  }
  if (is.list(value)) {
    if (!is.null(names(value))) {
      value <- value[order(names(value), method = "radix")]
    }
    return(lapply(value, .peptideQcCanonicalise))
  }
  value
}

.peptideQcDigest <- function(value) {
  digest::digest(
    .peptideQcCanonicalise(value),
    algo = "sha256",
    serialize = TRUE,
    serializeVersion = 2
  )
}

.peptideQcImplementationVersion <- function() {
  version <- tryCatch(
    as.character(utils::packageVersion("MultiScholaR")),
    error = function(e) NA_character_
  )
  if (is.na(version)) {
    description_file <- file.path(getwd(), "DESCRIPTION")
    version <- tryCatch(
      as.character(read.dcf(description_file, fields = "Version")[[1L]]),
      error = function(e) "development"
    )
  }
  paste0("MultiScholaR@", version)
}

.peptideQcObjectArgs <- function(object) {
  if (base::isS4(object) && "args" %in% methods::slotNames(object) &&
      is.list(object@args)) {
    return(object@args)
  }
  list()
}

.isPeptideQcAuditableObject <- function(object) {
  base::isS4(object) &&
    "args" %in% methods::slotNames(object) &&
    (methods::is(object, "PeptideQuantitativeData") ||
      methods::is(object, "ProteinQuantitativeData"))
}

.peptideQcObjectTable <- function(object) {
  if (!base::isS4(object)) {
    stop("Peptide-QC audit requires an S4 quantitative data object.", call. = FALSE)
  }
  if (methods::is(object, "PeptideQuantitativeData")) {
    return(object@peptide_data)
  }
  if (methods::is(object, "ProteinQuantitativeData")) {
    return(object@protein_quant_table)
  }
  stop(
    sprintf("Peptide-QC audit does not support S4 class `%s`.", class(object)[[1L]]),
    call. = FALSE
  )
}

.peptideQcColumnRoles <- function(object) {
  roles <- list(object_class = class(object)[[1L]])
  if (methods::is(object, "PeptideQuantitativeData")) {
    roles <- c(roles, list(
      active_protein_id = object@protein_id_column,
      stripped_peptide = object@peptide_sequence_column,
      modified_peptidoform = if ("Modified.Sequence" %in% names(object@peptide_data)) "Modified.Sequence" else NA_character_,
      run = object@sample_id,
      raw_quantity = object@raw_quantity_column,
      normalised_quantity = object@norm_quantity_column,
      q_value = object@q_value_column,
      global_q_value = object@global_q_value_column,
      global_pg_q_value = if ("Global.PG.Q.Value" %in% names(object@peptide_data)) "Global.PG.Q.Value" else NA_character_,
      proteotypic = object@proteotypic_peptide_sequence_column,
      biological_group = object@group_id,
      technical_replicate_group = object@technical_replicate_id
    ))
  } else {
    roles <- c(roles, list(
      active_protein_id = object@protein_id_column,
      run = object@sample_id,
      biological_group = object@group_id,
      technical_replicate_group = object@technical_replicate_id
    ))
  }
  roles
}

.peptideQcCountDistinct <- function(data, columns) {
  columns <- columns[!is.na(columns) & columns %in% names(data)]
  if (length(columns) == 0L || nrow(data) == 0L) return(0L)
  as.integer(nrow(unique(data[, columns, drop = FALSE])))
}

.peptideQcIdentificationSummary <- function(data, protein_column) {
  evidence_columns <- intersect(
    c("identification_peptide_count", "identification_peptidoform_count"),
    names(data)
  )
  if (length(evidence_columns) == 0L || !protein_column %in% names(data)) {
    return(list(status = "not_frozen", protein_count = 0L))
  }
  evidence <- unique(data[, c(protein_column, evidence_columns), drop = FALSE])
  evidence_range <- function(column) {
    if (!column %in% names(evidence) || nrow(evidence) == 0L) return(c(NA_real_, NA_real_))
    values <- suppressWarnings(as.numeric(evidence[[column]]))
    values <- values[!is.na(values) & is.finite(values)]
    if (length(values) == 0L) c(NA_real_, NA_real_) else range(values)
  }
  peptide_range <- evidence_range("identification_peptide_count")
  peptidoform_range <- evidence_range("identification_peptidoform_count")
  list(
    status = "frozen_after_identification_filter",
    protein_count = as.integer(nrow(evidence)),
    peptide_count_min = peptide_range[[1L]],
    peptide_count_max = peptide_range[[2L]],
    peptidoform_count_min = peptidoform_range[[1L]],
    peptidoform_count_max = peptidoform_range[[2L]]
  )
}

.peptideQcAuditSummary <- function(object) {
  data <- .peptideQcObjectTable(object)
  roles <- .peptideQcColumnRoles(object)
  protein_column <- roles$active_protein_id
  run_column <- roles$run
  stripped_column <- .peptideQcOr(roles$stripped_peptide, NA_character_)
  modified_column <- .peptideQcOr(roles$modified_peptidoform, NA_character_)

  quantity_columns <- unique(c(
    roles$raw_quantity,
    roles$normalised_quantity,
    "Peptide.Imputed",
    grep("(Quantity|Normalised|Imputed)$", names(data), value = TRUE)
  ))
  quantity_columns <- quantity_columns[
    !is.na(quantity_columns) & quantity_columns %in% names(data)
  ]
  quantity_columns <- quantity_columns[
    vapply(data[quantity_columns], is.numeric, logical(1))
  ]
  quantity_values <- if (length(quantity_columns) > 0L) {
    unlist(data[quantity_columns], use.names = FALSE)
  } else numeric(0)

  frozen_summary <- .peptideQcIdentificationSummary(data, protein_column)
  if (identical(frozen_summary$status, "not_frozen")) {
    prior_records <- .peptideQcObjectArgs(object)$peptide_qc_audit$records
    prior_frozen <- lapply(.peptideQcOr(prior_records, list()), function(record) {
      record$after_summary$frozen_identification_evidence
    })
    prior_frozen <- Filter(
      function(summary) is.list(summary) && identical(summary$status, "frozen_after_identification_filter"),
      prior_frozen
    )
    if (length(prior_frozen) > 0L) frozen_summary <- tail(prior_frozen, 1L)[[1L]]
  }

  current_survivors <- list(
    active_protein_groups = .peptideQcCountDistinct(data, protein_column),
    protein_stripped_peptide_pairs = .peptideQcCountDistinct(data, c(protein_column, stripped_column)),
    protein_modified_peptidoform_pairs = .peptideQcCountDistinct(
      data,
      c(protein_column, if (!is.na(modified_column)) modified_column else stripped_column)
    )
  )
  args <- .peptideQcObjectArgs(object)
  exclusion_summary <- args$srlQvalueProteotypicPeptideClean$biological_exclusion_summary
  excluded <- if (is.data.frame(exclusion_summary) && nrow(exclusion_summary) > 0L &&
      "excluded_rows" %in% names(exclusion_summary)) {
    as.integer(exclusion_summary$excluded_rows[[1L]])
  } else 0L

  list(
    rows = as.integer(nrow(data)),
    distinct_runs = .peptideQcCountDistinct(data, run_column),
    active_protein_groups = current_survivors$active_protein_groups,
    distinct_stripped_peptides = .peptideQcCountDistinct(data, stripped_column),
    distinct_modified_peptidoforms = .peptideQcCountDistinct(
      data,
      if (!is.na(modified_column)) modified_column else stripped_column
    ),
    frozen_identification_evidence = frozen_summary,
    current_quantitative_survivors = current_survivors,
    quantity_columns = quantity_columns,
    missing_quantity_values = as.integer(sum(is.na(quantity_values))),
    non_finite_quantity_values = as.integer(sum(!is.na(quantity_values) & !is.finite(quantity_values))),
    imputed_cells = if ("is_imputed" %in% names(data)) as.integer(sum(data$is_imputed %in% TRUE, na.rm = TRUE)) else 0L,
    decoy_rows = if ("is_decoy" %in% names(data)) as.integer(sum(data$is_decoy %in% TRUE, na.rm = TRUE)) else 0L,
    contaminant_rows = if ("is_contaminant" %in% names(data)) as.integer(sum(data$is_contaminant %in% TRUE, na.rm = TRUE)) else 0L,
    biologically_excluded_rows = excluded
  )
}

.peptideQcStableIdentities <- function(object) {
  data <- .peptideQcObjectTable(object)
  roles <- .peptideQcColumnRoles(object)
  candidates <- unique(c(
    roles$active_protein_id,
    roles$run,
    roles$stripped_peptide,
    roles$modified_peptidoform,
    "Precursor.Id",
    "Precursor.Charge"
  ))
  candidates <- candidates[!is.na(candidates) & candidates %in% names(data)]
  identities <- data[, candidates, drop = FALSE]
  identities[] <- lapply(identities, function(value) {
    if (is.factor(value)) as.character(value) else value
  })
  unique(identities)
}

.peptideQcIdentityKey <- function(data, columns = names(data)) {
  if (nrow(data) == 0L) return(character(0))
  columns <- columns[columns %in% names(data)]
  do.call(paste, c(lapply(data[columns], function(value) {
    value <- ifelse(is.na(value), "<NA>", as.character(value))
    paste0(nchar(value, type = "bytes"), ":", value)
  }), sep = "\u001f"))
}

.peptideQcConfidenceFailureReasons <- function(data, params) {
  q_columns <- c("Q.Value", "Global.Q.Value", "Global.PG.Q.Value")
  thresholds <- c(
    .peptideQcOr(params$qvalue_threshold, 0.01),
    .peptideQcOr(params$global_qvalue_threshold, 0.01),
    .peptideQcOr(params$global_pg_qvalue_threshold, 0.01)
  )
  reasons <- rep("", nrow(data))
  for (index in seq_along(q_columns)) {
    column <- q_columns[[index]]
    if (!column %in% names(data)) next
    values <- suppressWarnings(as.numeric(data[[column]]))
    missing <- is.na(values) | !is.finite(values)
    above <- !missing & values > thresholds[[index]]
    reasons[missing] <- paste0(reasons[missing], ifelse(nzchar(reasons[missing]), ";", ""), column, "_missing_or_non_finite")
    threshold_label <- format(thresholds[[index]], scientific = FALSE, trim = TRUE)
    reasons[above] <- paste0(
      reasons[above],
      ifelse(nzchar(reasons[above]), ";", ""),
      column,
      "_above_",
      threshold_label
    )
  }
  if (isTRUE(params$proteotypic_only) && "Proteotypic" %in% names(data)) {
    failed <- is.na(data$Proteotypic) | suppressWarnings(as.numeric(data$Proteotypic)) != 1
    reasons[failed] <- paste0(reasons[failed], ifelse(nzchar(reasons[failed]), ";", ""), "not_proteotypic")
  }
  reasons[!nzchar(reasons)] <- "excluded_by_qvalue_or_biological_classification"
  reasons
}

.peptideQcSpecialisedRemovalLedger <- function(object, stage_id) {
  args <- .peptideQcObjectArgs(object)
  section <- switch(
    stage_id,
    intensity_filter = args$peptideIntensityFiltering,
    protein_evidence_filter = args$filterMinNumPeptidesPerProtein,
    sample_filter = args$filterMinNumPeptidesPerSample,
    replicate_filter = args$removePeptidesWithOnlyOneReplicate,
    qvalue_filter = args$srlQvalueProteotypicPeptideClean,
    NULL
  )
  ledger <- if (identical(stage_id, "qvalue_filter")) {
    section$biological_exclusion_ledger
  } else {
    section$removal_ledger
  }
  if (is.data.frame(ledger)) ledger else NULL
}

.peptideQcRemovalLedger <- function(before, after, stage_id, params) {
  before_ids <- .peptideQcStableIdentities(before)
  after_ids <- .peptideQcStableIdentities(after)
  shared <- intersect(names(before_ids), names(after_ids))
  if (length(shared) == 0L || nrow(before_ids) == 0L) {
    return(data.frame())
  }
  after_keys <- .peptideQcIdentityKey(after_ids, shared)
  removed <- before_ids[!.peptideQcIdentityKey(before_ids, shared) %in% after_keys, , drop = FALSE]
  if (nrow(removed) == 0L) return(removed)

  if (identical(stage_id, "qvalue_filter")) {
    source_data <- .peptideQcObjectTable(before)
    source_key <- .peptideQcIdentityKey(source_data, intersect(names(source_data), names(removed)))
    reason_by_key <- split(
      .peptideQcConfidenceFailureReasons(source_data, params),
      source_key
    )
    removed_key <- .peptideQcIdentityKey(removed, intersect(names(source_data), names(removed)))
    removed$failure_reason <- vapply(removed_key, function(key) {
      paste(sort(unique(.peptideQcOr(reason_by_key[[key]], "excluded_by_qvalue_or_biological_classification"))), collapse = ";")
    }, character(1))
  } else if (identical(stage_id, "protein_evidence_filter")) {
    before_data <- .peptideQcObjectTable(before)
    protein_column <- .peptideQcColumnRoles(before)$active_protein_id
    evidence_columns <- c(
      protein_column,
      "identification_peptide_count",
      "identification_peptidoform_count"
    )
    if (all(evidence_columns %in% names(before_data))) {
      evidence <- unique(before_data[, evidence_columns, drop = FALSE])
      removed <- merge(removed, evidence, by = protein_column, all.x = TRUE, sort = FALSE)
    }
    peptide_threshold <- as.numeric(
      .peptideQcOr(
        params$num_peptides_per_protein_thresh,
        .peptideQcOr(
          params$min_peptides_per_protein,
          .peptideQcOr(params$peptides_per_protein_cutoff, 1)
        )
      )
    )
    peptidoform_threshold <- as.numeric(
      .peptideQcOr(
        params$num_peptidoforms_per_protein_thresh,
        .peptideQcOr(
          params$min_peptidoforms_per_protein,
          .peptideQcOr(params$peptidoforms_per_protein_cutoff, 2)
        )
      )
    )
    peptide_failed <- if ("identification_peptide_count" %in% names(removed)) {
      is.na(removed$identification_peptide_count) |
        removed$identification_peptide_count < peptide_threshold
    } else rep(TRUE, nrow(removed))
    peptidoform_failed <- if ("identification_peptidoform_count" %in% names(removed)) {
      is.na(removed$identification_peptidoform_count) |
        removed$identification_peptidoform_count < peptidoform_threshold
    } else rep(TRUE, nrow(removed))
    removed$required_identification_peptide_count <- peptide_threshold
    removed$required_identification_peptidoform_count <- peptidoform_threshold
    removed$failure_reason <- ifelse(
      peptide_failed & peptidoform_failed,
      "below_peptide_and_peptidoform_thresholds",
      ifelse(peptide_failed, "below_peptide_threshold", "below_peptidoform_threshold")
    )
  } else {
    removed$failure_reason <- switch(
      stage_id,
      intensity_filter = "failed_intensity_and_group_support_contract",
      protein_evidence_filter = "below_frozen_identification_evidence_threshold",
      sample_filter = "sample_below_distinct_peptide_identity_threshold",
      replicate_filter = "no_replicate_group_with_two_distinct_runs",
      "removed_by_stage_contract"
    )
  }

  specialised <- .peptideQcSpecialisedRemovalLedger(after, stage_id)
  if (is.data.frame(specialised) && nrow(specialised) > 0L) {
    join_columns <- intersect(names(removed), names(specialised))
    join_columns <- setdiff(join_columns, "failure_reason")
    if (length(join_columns) > 0L) {
      names(specialised)[names(specialised) == "failure_reason"] <- "stage_failure_reason"
      removed <- merge(removed, specialised, by = join_columns, all.x = TRUE, sort = FALSE)
      if ("stage_failure_reason" %in% names(removed)) {
        use_stage <- !is.na(removed$stage_failure_reason) & nzchar(removed$stage_failure_reason)
        removed$failure_reason[use_stage] <- removed$stage_failure_reason[use_stage]
      }
      if ("exclusion_reason" %in% names(removed)) {
        use_exclusion <- !is.na(removed$exclusion_reason) & nzchar(removed$exclusion_reason)
        removed$failure_reason[use_exclusion] <- removed$exclusion_reason[use_exclusion]
      }
    }
  }
  removed$source_stage <- stage_id
  removed
}

.peptideQcImputationLedger <- function(object, stage_id) {
  data <- .peptideQcObjectTable(object)
  if (!"is_imputed" %in% names(data) || !any(data$is_imputed %in% TRUE, na.rm = TRUE)) {
    return(data.frame())
  }
  roles <- .peptideQcColumnRoles(object)
  identity_columns <- intersect(
    unique(c(roles$active_protein_id, roles$run, roles$stripped_peptide, roles$modified_peptidoform)),
    names(data)
  )
  value_columns <- intersect(c("Peptide.Imputed", roles$normalised_quantity), names(data))
  ledger <- data[data$is_imputed %in% TRUE, unique(c(identity_columns, value_columns)), drop = FALSE]
  ledger$source_stage <- stage_id
  ledger$imputation_reason <- "eligible_missing_value_with_observed_technical_replicate_support"
  unique(ledger)
}

.validatePeptideQcStagePrerequisites <- function(object, stage_id) {
  data <- .peptideQcObjectTable(object)
  roles <- .peptideQcColumnRoles(object)
  if (identical(stage_id, "protein_evidence_filter")) {
    frozen <- c("identification_peptide_count", "identification_peptidoform_count")
    legacy_identity <- c(roles$stripped_peptide, "Modified.Sequence")
    legacy_rollup <- c(roles$stripped_peptide, "peptidoform_ids")
    has_evidence <- all(frozen %in% names(data))
    has_legacy_identity <- all(legacy_identity %in% names(data))
    has_legacy_rollup <- all(legacy_rollup %in% names(data))
    if (!has_evidence && !has_legacy_identity && !has_legacy_rollup) {
      stop(
        paste0(
          "Peptide-QC stage `protein_evidence_filter` prerequisite error: ",
          "requires frozen identification counts or both stripped and modified ",
          "peptide identities."
        ),
        call. = FALSE
      )
    }
    return(list(
      required_columns = c(roles$active_protein_id, roles$stripped_peptide),
      evidence_contract = if (has_evidence) "frozen_identification_counts" else "legacy_identity_fallback",
      validated = TRUE
    ))
  }
  required <- switch(
    stage_id,
    qvalue_filter = c(roles$active_protein_id, roles$run, roles$stripped_peptide, "Modified.Sequence", "Q.Value", "Global.Q.Value", "Global.PG.Q.Value", "Proteotypic"),
    precursor_rollup = c(roles$active_protein_id, roles$run, roles$stripped_peptide, "Modified.Sequence"),
    intensity_filter = c(roles$active_protein_id, roles$run, roles$stripped_peptide),
    sample_filter = c(roles$active_protein_id, roles$run, roles$stripped_peptide),
    replicate_filter = c(roles$active_protein_id, roles$run, roles$stripped_peptide),
    imputation = c(roles$active_protein_id, roles$run, roles$stripped_peptide),
    protein_rollup = roles$active_protein_id,
    character(0)
  )
  required <- unique(required[!is.na(required)])
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop(
      sprintf(
        "Peptide-QC stage `%s` prerequisite error: missing required column(s): %s.",
        stage_id,
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  list(required_columns = required, validated = TRUE)
}

.peptideQcCurrentStateName <- function(state_manager) {
  current <- workflowStateCurrentName(state_manager)
  if (!is.null(current) && length(current) == 1L && nzchar(current)) return(current)
  history <- tryCatch(state_manager$getHistory(), error = function(e) character(0))
  if (length(history) > 0L) tail(history, 1L) else "legacy_untracked_state"
}

.peptideQcImmutableImport <- function(object) {
  args <- .peptideQcObjectArgs(object)
  existing <- args$peptide_qc_audit$immutable_import
  if (is.list(existing) && !is.null(existing$digest)) return(existing)
  list(
    schema_version = .peptideQcAuditSchemaVersion,
    digest = .peptideQcDigest(list(
      data = .peptideQcObjectTable(object),
      column_roles = .peptideQcColumnRoles(object)
    )),
    source_label = "in_memory_import",
    source_path_recorded = FALSE
  )
}

.emitPeptideQcAuditRecord <- function(before,
                                       after,
                                       stage_id,
                                       parent_state,
                                       current_state,
                                       resolved_parameters = list(),
                                       status = "applied",
                                       decision_reason = NA_character_,
                                       now = Sys.time(),
                                       transformation_type = "filter") {
  if (!base::isS4(after) || !"args" %in% methods::slotNames(after)) return(after)
  if (!status %in% c("applied", "skipped", "no_op")) {
    stop("Audit status must be applied, skipped, or no_op.", call. = FALSE)
  }

  before_args <- .peptideQcObjectArgs(before)
  lineage <- before_args$peptide_qc_audit
  immutable_import <- .peptideQcImmutableImport(before)
  prior_records <- lineage$records
  if (!is.list(prior_records)) prior_records <- list()
  parent_record_id <- if (length(prior_records) > 0L) {
    prior_records[[length(prior_records)]]$record_id
  } else {
    paste0("legacy-import:", substr(immutable_import$digest, 1L, 16L))
  }

  prereqs <- if (identical(status, "skipped")) {
    list(validated = FALSE, skip_reason = decision_reason)
  } else {
    .validatePeptideQcStagePrerequisites(before, stage_id)
  }
  before_summary <- .peptideQcAuditSummary(before)
  after_summary <- .peptideQcAuditSummary(after)
  removal_ledger <- if (identical(transformation_type, "filter")) {
    .peptideQcRemovalLedger(before, after, stage_id, resolved_parameters)
  } else data.frame()
  imputation_ledger <- .peptideQcImputationLedger(after, stage_id)

  substantive <- list(
    schema_version = .peptideQcAuditSchemaVersion,
    stage_id = stage_id,
    status = status,
    decision_reason = decision_reason,
    transformation_type = transformation_type,
    parent_state = parent_state,
    current_state = current_state,
    parent_record_id = parent_record_id,
    implementation_version = .peptideQcImplementationVersion(),
    resolved_parameters = resolved_parameters,
    column_roles = .peptideQcColumnRoles(after),
    prerequisites = prereqs,
    before_summary = before_summary,
    after_summary = after_summary,
    removal_ledger = removal_ledger,
    imputation_ledger = imputation_ledger,
    immutable_import_digest = immutable_import$digest,
    data_digest = .peptideQcDigest(.peptideQcObjectTable(after)),
    configuration_digest = .peptideQcDigest(resolved_parameters)
  )
  record_id <- paste0("pqc:", substr(.peptideQcDigest(substantive), 1L, 24L))
  record <- c(
    substantive,
    list(
      record_id = record_id,
      timestamp_utc = format(now, "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC"),
      canonical_digest = .peptideQcDigest(substantive)
    )
  )
  prior_records[[length(prior_records) + 1L]] <- record
  after@args$peptide_qc_audit <- list(
    schema_version = .peptideQcAuditSchemaVersion,
    provenance_status = if (length(prior_records) == 1L) "legacy_input_wrapped" else "tracked",
    immutable_import = immutable_import,
    records = prior_records,
    current_record_id = record_id
  )
  after
}

.savePeptideQcState <- function(state_manager,
                                 before,
                                 after,
                                 stage_id,
                                 state_name,
                                 config_object,
                                 description,
                                 audit_parameters = config_object,
                                 now = Sys.time(),
                                 status = "applied",
                                 decision_reason = NA_character_,
                                 transformation_type = "filter") {
  audit_enabled <- workflowStateAuditEnabled(state_manager)
  if (!identical(audit_enabled, FALSE) &&
      .isPeptideQcAuditableObject(before) &&
      .isPeptideQcAuditableObject(after)) {
    after <- .emitPeptideQcAuditRecord(
      before = before,
      after = after,
      stage_id = stage_id,
      parent_state = .peptideQcCurrentStateName(state_manager),
      current_state = state_name,
      resolved_parameters = audit_parameters,
      status = status,
      decision_reason = decision_reason,
      now = now,
      transformation_type = transformation_type
    )
    config_object$audit_record_id <- after@args$peptide_qc_audit$current_record_id
  }
  state_manager$saveState(
    state_name = state_name,
    s4_data_object = after,
    config_object = config_object,
    description = description
  )
  after
}

.peptideQcSubstantiveAudit <- function(object) {
  audit <- .peptideQcObjectArgs(object)$peptide_qc_audit
  if (!is.list(audit)) return(list(provenance_status = "legacy_untracked"))
  records <- lapply(audit$records, function(record) {
    record$timestamp_utc <- NULL
    record
  })
  list(
    schema_version = audit$schema_version,
    provenance_status = audit$provenance_status,
    immutable_import = audit$immutable_import,
    records = records,
    current_record_id = audit$current_record_id
  )
}
