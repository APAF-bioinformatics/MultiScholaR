# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# KIT-108 W4 (D3): proteomics' route pieces for the shared Summary & Export tab (mod_omic_summary.R, the core in
# utils_workflow_summary.R), named by `omicGuiConfig("proteomics")$summary`. The copy's gate, the integration file's
# name and its in-process writer are today's (mod_prot_summary_publication_helpers.R:104;
# mod_prot_summary_session_helpers.R:82-98 at 1cb2a79b), and the session-state payload's adapter calls today's
# `prepareProtSummarySessionStateExport` with the core's arguments. Below them, the seven functions D3's map moves here
# unchanged by name and body (mod_prot_summary_report_helpers.R:1-248, mod_prot_summary_session_helpers.R:144-215; the
# files are deleted by W4): the final state, the four-template status, the route's template resolver family and the
# session-state payload. F104-01 (D5): the native export's metadata (`collectDIANativeSessionMetadata`) and the
# summary's archive refresh (`refreshDIANativeSummaryArchive`) name the template `resolveProtSummaryReportTemplate`
# resolves, as the report step does. KIT-061: the integration file is the integration input's folder
# (utils_workflow_integration_input.R), written from the route's final matrix; the native archive and the S4 RDS are no
# longer written to the integration folder (KIT-060 D2).

#' The copy's gate: the session's owner has no export pending (today's `assertDIAExportReady(getDIAActiveManager())`).
#' @noRd
protSummaryCopyGate <- function() assertDIAExportReady(getDIAActiveManager())

#' The integration input's manifest (KIT-061, D1 and D2): `<omic>/integration_input.json` under the project's
#' integration folder, a fixed name with no label.
#' @noRd
protSummaryIntegrationFilename <- function(workflowData, omicType, label) {
  file.path(omicType, INTEGRATION_INPUT_MANIFEST)
}

#' The integration input, written in process (KIT-061): the contract's folder (`writeIntegrationInput`) from the
#' route's final matrix (`protIntegrationInputSource`). On a native session the summary's archive in the source
#' directory is still refreshed (the session export, which the publication copy carries), but no longer copied to the
#' integration folder (D2); no S4 object is written.
#' @param path The manifest's path (`protSummaryIntegrationFilename`'s); its folder is the omic's input.
#' @param object,saveRDSFn Unused: the route writer's arguments in the core's call.
#' @param sourceDir The project's source directory (the native archive's and the UniProt annotation's).
#' @return The notice the save shows, or FALSE when the owner holds no final matrix (K061-3: an earlier save's input
#'   is then removed).
#' @noRd
protSummaryIntegrationWrite <- function(workflowData, path, object = NULL, sourceDir = NULL, saveRDSFn = saveRDS) {
  manager <- workflowData$state_manager
  written <- FALSE
  if (inherits(manager, "DIAStateAdapter")) {
    if (isDIANativeSession(manager)) refreshDIANativeSummaryArchive(workflowData, sourceDir)
    annotation <- workflowData$uniprot_dat_cln
    annotationFile <- if (!is.null(sourceDir)) file.path(sourceDir, "uniprot_dat_cln.RDS") else ""
    if (is.null(annotation) && file.exists(annotationFile)) annotation <- readRDS(annotationFile)
    written <- withDIAStateScope(manager, function() {
      source <- protIntegrationInputSource(manager)
      if (is.null(source)) return(FALSE)
      writeIntegrationInput(dirname(path), source$manifest, source$samples,
        matrices = list(proteomics = source$matrix),
        features = list(proteomics = protIntegrationInputFeatures(source$matrix$feature_id, annotation,
          mapping = if (referenceMappingEnabled()) readReferenceMapping(manager)$mapping)),
        corrected = if (!is.null(source$corrected)) list(proteomics = source$corrected) else list())
      TRUE
    })
  }
  if (isTRUE(written)) return("Saved Integration input")
  removeIntegrationInput(dirname(path))
  FALSE
}

#' The head's lineage (K061-1): the current state, then each record's `previous_state`, until a state without a
#' reference, a repeated name, a state stored as a whole legacy object (its record NULL) or one of `until`.
#' @param manager The session's owner, inside its state scope.
#' @param until The state names that end the walk (included).
#' @return Named list of the states' records, head first.
#' @noRd
protIntegrationInputLineage <- function(manager, until) {
  lineage <- list()
  state <- manager$current_state
  while (is.character(state) && length(state) == 1L && !is.na(state) && !state %in% names(lineage)) {
    reference <- getDIAStateReference(manager, state)
    if (is.null(reference)) break
    # Every declared class, as the analysis record reads the lineage (KIT-101 F-15): a result-class head (the
    # proteomics `enrichment_completed`) is walked through.
    record <- tryCatch(readComponentMetadata(manager, reference), saved_entity_unsupported = function(error) NULL)
    lineage[state] <- list(record)
    if (is.null(record) || state %in% until) break
    state <- record$previous_state
  }
  lineage
}

#' The proteomics integration input's matrix and manifest fields (KIT-061; KIT-060 audit 4.3, D14): the route's final
#' state nearest the head in the head's own lineage (K061-1: the current state when it is one; `normalized` or
#' `ruv_corrected`, then `correlation_filtered` on DIA, LFQ and TMT; `correlation_filtered`, then `ruv_covariates` or
#' `ruv_corrected` on DIA_LIMPA), never a name kept in the history from an earlier run. Its matrix is the state's
#' protein table, except on `ruv_covariates`, whose table is uncorrected: there the RUV-adjusted product. The role is
#' `corrected` when the lineage's normalisation or RUV state is an RUV state with RUV applied, else `normalised`; on
#' DIA_LIMPA a `correlation_filtered` final is always `normalised` (its input is the pre-RUV matrix). KIT-153 (INT-23):
#' where RUV was applied, that matrix is the `corrected` one, and `matrix` is the normalised matrix the RUV step read
#' (`protIntegrationInputRuvInput`) on its proteins and samples; the manifest's role is always `normalised`.
#' @param manager The session's owner, inside its state scope.
#' @return `list(manifest, samples, matrix, corrected)` (each matrix's first column `feature_id`, its samples in the
#'   design's order; `corrected` NULL without RUV), or NULL when the lineage has no final state or it is stored as a
#'   whole legacy object.
#' @noRd
protIntegrationInputSource <- function(manager) {
  ruvStates <- c("ruv_corrected", "ruv_covariates")
  # The walk ends at the first normalisation or RUV state: it decides the role.
  lineage <- protIntegrationInputLineage(manager, until = c("normalized", ruvStates))
  final <- utils::head(intersect(names(lineage), c("normalized", "ruv_corrected", "correlation_filtered", "ruv_covariates")), 1L)
  if (!length(final)) return(NULL)
  reference <- getDIAStateReference(manager, final)
  # A state stored as a whole legacy object has no component record; it is not hydrated here (as the analysis
  # record leaves it), so no input is written for it.
  record <- lineage[[final]]$data
  if (is.null(record)) return(NULL)
  workflow <- manager$workflow_type %||% manager$.context$snapshot$workflow_type %||% NA_character_
  decider <- utils::tail(names(lineage), 1L)
  ruvState <- if (decider %in% ruvStates && !identical(lineage[[decider]]$config$ruv_applied, FALSE)) decider else character()
  # A record's previous state is the head at its commit, not its input, and no checkpoint records its source. On
  # DIA_LIMPA the correlation filter always reads the pre-RUV input (`diaCorrelationSource`), and the route corrects
  # only in ruv_covariates' adjusted product or ruv_corrected, so its correlation_filtered table is normalised. The
  # classic routes' correlation reads the RUV output, so there it stays corrected after RUV.
  if (identical(workflow, "DIA_LIMPA") && identical(final, "correlation_filtered")) ruvState <- character()
  slots <- record$slots
  matrixSource <- reference
  product <- NULL
  # KIT-129 (audit F2): a `ruv_covariates` state committed without RUV has no adjusted product; its own table is
  # the matrix.
  if (identical(final, "ruv_covariates") && length(ruvState)) {
    science <- hydrateSavedWorkflowEntity(manager$.context$store, record$products$ruv_covariates_science)$value
    matrixSource <- science$adjusted
    product <- list(type = "ruv-adjusted-science", entity_id = matrixSource$entity_id,
      sha256 = matrixSource$payload$sha256)
  }
  table <- readDIAStateMetadata(manager, matrixSource)$data$components$protein_quant_table
  # The class's validity needs the design beside the table.
  matrix <- readDIAExecutionComponents(manager, matrixSource,
    fields = c("protein_quant_table", "design_matrix"))@protein_quant_table
  design <- readDIAExecutionComponents(manager, reference, fields = "design_matrix")@design_matrix
  samples <- setdiff(names(matrix), slots$protein_id_column)
  design <- design[match(samples, as.character(design[[slots$sample_id]])), , drop = FALSE]
  rownames(design) <- NULL
  matrix <- data.frame(feature_id = as.character(matrix[[slots$protein_id_column]]), matrix[samples],
    check.names = FALSE, stringsAsFactors = FALSE)
  # KIT-153 (INT-23): after RUV the matrix above is the corrected one; the matrix written as the data set's own is the
  # normalised matrix the RUV step read (`protIntegrationInputRuvInput`), on the corrected matrix's proteins and samples.
  corrected <- NULL
  normalised <- list(name = final, entity_id = reference$entity_id, sha256 = reference$payload$sha256)
  if (length(ruvState)) {
    input <- protIntegrationInputRuvInput(manager, lineage[[ruvState]])
    uncorrected <- readDIAExecutionComponents(manager, input$reference,
      fields = c("protein_quant_table", "design_matrix"))@protein_quant_table
    uncorrected <- data.frame(feature_id = as.character(uncorrected[[slots$protein_id_column]]),
      uncorrected[setdiff(names(uncorrected), slots$protein_id_column)], check.names = FALSE, stringsAsFactors = FALSE)
    corrected <- matrix
    matrix <- integrationInputUncorrectedOn(uncorrected, corrected, "proteomics")
    uncorrected <- NULL
    normalised <- list(name = input$name, entity_id = input$reference$entity_id, sha256 = input$reference$payload$sha256)
  }
  manifest <- list(omic_type = "proteomics", workflow = workflow,
    matrix_role = "normalised", scale = "log2",
    normalisation = list(method = protIntegrationInputNormMethod(manager), source = normalised),
    correction = if (length(ruvState)) protIntegrationInputCorrection(manager, ruvState, design) else NULL,
    source_state = list(name = final, entity_id = reference$entity_id, sha256 = reference$payload$sha256,
      protein_id_column = slots$protein_id_column,
      matrix = list(entity_id = table$entity_id, sha256 = table$payload$sha256), product = product),
    sample_id_column = slots$sample_id, group_column = slots$group_id,
    technical_replicate_column = slots$technical_replicate_id %||% NA_character_,
    # KIT-062: the biological sample's column, by the rule the other omics' writers use (INT-2).
    biological_sample_column = integrationInputBiologicalSampleColumn(design,
      slots$technical_replicate_id %||% NA_character_, slots$sample_id))
  list(manifest = manifest, samples = design, matrix = matrix, corrected = corrected)
}

#' The normalised matrix a proteomics RUV step read (KIT-153, INT-23), from the RUV state's own product: the classic
#' route's normalisation and RUV-III-C (`normalization_science`, its `normalized` checkpoint: the between-sample
#' normalisation the RUV kernel corrected), or DIA_LIMPA's RUV-III fit (`ruv_covariates_science`, its `covariates`
#' checkpoint: the step's input table, with W1..Wk in the design only). Refused in words when the product does not
#' record it.
#' @param manager The session's owner, inside its state scope.
#' @param record The RUV state's record (the lineage's).
#' @return `list(name, reference)`.
#' @noRd
protIntegrationInputRuvInput <- function(manager, record) {
  products <- record$data$products
  covariates <- !is.null(products$ruv_covariates_science)
  product <- if (covariates) products$ruv_covariates_science else products$normalization_science
  # As `protIntegrationInputCorrection` reads it: a product hydrates as an entity whose `value` is the kernel's result; a
  # plain saved value is the list itself.
  hydrated <- if (!is.null(product)) hydrateSavedWorkflowEntity(manager$.context$store, product)
  value <- if (is.list(hydrated$value)) hydrated$value else hydrated
  input <- if (is.list(value)) (if (covariates) value$covariates else value$normalized)
  if (!isSavedWorkflowReference(input)) {
    abortIntegrationInput(paste("the normalised matrix that entered RUV is not recorded with the RUV step; run the",
      "normalisation again, then save"))
  }
  list(name = if (covariates) "ruv_step_input" else "normalized", reference = input)
}

#' The between-sample normalisation of a proteomics input (KIT-153): the `norm_method` of the state nearest the head
#' in its lineage that records one, else NA.
#' @param manager The session's owner, inside its state scope.
#' @noRd
protIntegrationInputNormMethod <- function(manager) {
  for (record in protIntegrationInputLineage(manager, until = character())) {
    method <- record$config$norm_method
    if (is.character(method) && length(method) == 1L && !is.na(method)) return(method)
  }
  NA_character_
}

#' The RUV design that produced a corrected proteomics matrix: the classic route's RUV-III-C (KIT-060 audit row C1:
#' the normalisation's commit and product), or DIA_LIMPA's RUV-III fit (row C2: the covariates step's commit, in either
#' mode).
#' @param state The RUV state's name.
#' @param design The samples' design (for the grouping variable's levels).
#' @return A list for the manifest's `correction`.
#' @noRd
protIntegrationInputCorrection <- function(manager, state, design) {
  record <- readDIAStateMetadata(manager, getDIAStateReference(manager, state))
  config <- record$config
  products <- record$data$products
  science <- function(product) {
    # As the record reads it: a product hydrates as an entity whose `value` is the kernel's result; a plain saved value
    # is the list itself.
    hydrated <- hydrateSavedWorkflowEntity(manager$.context$store, product)
    if (is.list(hydrated$value)) hydrated$value else hydrated
  }
  correction <- if (!is.null(products$ruv_covariates_science)) {
    list(id = "C2", method = "RUVIII", ruv_mode = config$ruv_mode, k = config$ruv_k,
      percentage_as_neg_ctrl = config$percentage_as_neg_ctrl, control_count = config$control_count,
      replicate_structure = config$replicate_structure, grouping_variable = config$ruv_grouping_variable)
  } else {
    value <- if (!is.null(products$normalization_science)) science(products$normalization_science) else list()
    list(id = "C1", method = "RUVIII_C_Varying", ruv_mode = config$ruv_mode, k = config$ruv_k,
      percentage_as_neg_ctrl = config$percentage_as_neg_ctrl,
      control_count = if (is.logical(value$controls)) sum(value$controls, na.rm = TRUE) else NA_integer_,
      replicate_structure = record$data$slots$args$ruvIII_C_Varying$ruv_replicate_structure,
      grouping_variable = value$parameters$grouping)
  }
  # KIT-129 seventh pass item 3: a new design's run also names its design, the k applied, the replicate structure and
  # the rest of the step record's fields (`RUV_DESIGN_RECORD_FIELDS`, as the analysis record reads them). A legacy run
  # records none, and its correction is as before.
  operation <- if (!is.null(products$ruv_covariates_science)) "protein_ruv_covariates" else "normalization-ruv"
  outcomes <- workflowRecordStepOutcomes(manager$.context, state, operation, record, list())
  if (all(RUV_DESIGN_RECORD_FIELDS %in% names(outcomes))) correction <- c(correction, outcomes[RUV_DESIGN_RECORD_FIELDS])
  correction <- lapply(correction, function(value) if (is.null(value)) NA else value)
  grouping <- correction$grouping_variable
  correction$grouping_levels <- I(if (is.character(grouping) && grouping %in% names(design))
    unique(as.character(design[[grouping]])) else character())
  c(list(state = state), correction)
}

#' The proteomics feature annotation (D11: the omic's own): the UniProt accession (K061-4: the group's first member,
#' without a trailing isoform suffix `-<digits>`) and its `gene_names`; the display name is the first gene name, else
#' the feature id, which stays the key.
#' @param ids The matrix's feature ids.
#' @param annotation The UniProt table (`uniprot_dat_cln`, keyed by `Entry`), or NULL.
#' @param mapping NULL, or the project's reference mapping table (KIT-167): the features then gain
#'   `reference_accession_used` (by the protein-group rule, 8.3 items 1 and 8), `reference_gene`, `locus_tag` (NA: the
#'   mapping holds no strain locus tag) and `map_status` (`mapped`, `ambiguous_group` or `unmapped`); without a mapping
#'   the table is unchanged.
#' @return `data.frame(feature_id, display_name, uniprot_accession, gene_names)`, and the four mapping columns.
#' @noRd
protIntegrationInputFeatures <- function(ids, annotation, mapping = NULL) {
  accession <- sub("-[0-9]+$", "", sub(";.*$", "", ids))
  genes <- rep(NA_character_, length(ids))
  if (is.data.frame(annotation) && all(c("Entry", "gene_names") %in% names(annotation))) {
    genes <- as.character(annotation$gene_names)[match(accession, as.character(annotation$Entry))]
  }
  first <- vapply(strsplit(ifelse(is.na(genes), "", genes), " |:"), function(parts) {
    if (length(parts)) parts[[1L]] else ""
  }, "")
  features <- data.frame(feature_id = ids, display_name = ifelse(nzchar(first), first, ids), uniprot_accession = accession,
    gene_names = genes, stringsAsFactors = FALSE)
  if (is.null(mapping)) return(features)
  groups <- referenceAccessionForGroups(ids, mapping)
  features$reference_accession_used <- groups$reference_accession_used
  gene <- if ("reference_gene" %in% names(mapping)) mapping$reference_gene else rep(NA_character_, nrow(mapping))
  features$reference_gene <- gene[match(groups$reference_accession_used, mapping$reference_accession_used)]
  features$locus_tag <- rep(NA_character_, length(ids))
  features$map_status <- groups$group_class
  features
}

#' The session-state payload in the core's shape (`session_state`): today's `prepareProtSummarySessionStateExport`'s
#' state (the core writes it at the same path).
#' @noRd
protSummarySessionState <- function(inputValues, projectDirs, omicType, values, timestamp, workflowData = NULL) {
  prepareProtSummarySessionStateExport(projectDirs = projectDirs, omicType = omicType,
    experimentLabel = inputValues$experiment_label, description = inputValues$description,
    workflowArgsSaved = values$workflow_args_saved, filesCopied = values$files_copied,
    reportGenerated = values$report_generated, reportPath = values$report_path, workflowData = workflowData,
    timestamp = timestamp)$sessionState
}

# ------------------------------------------------------------------------------------------------
# Moved unchanged (D3's map): mod_prot_summary_report_helpers.R:1-248 and mod_prot_summary_session_helpers.R:144-215.
# ------------------------------------------------------------------------------------------------

resolveProtSummaryFinalS4State <- function(workflowData,
                                           logPrefix = "SESSION SUMMARY",
                                           catFn = cat) {
  finalS4Object <- NULL
  dataStateUsed <- NULL
  stateManager <- workflowData$state_manager

  if (is.null(stateManager)) {
    catFn(sprintf("%s: No state manager available\n", logPrefix))
    return(list(finalS4Object = NULL, dataStateUsed = NULL))
  }

  dataStates <- c(
    "ruv_covariates",
    "correlation_filtered",
    "ruv_corrected",
    "protein_replicate_filtered",
    "imputed"
  )
  availableStates <- diaStateHistory(stateManager)
  dataStateUsed <- purrr::detect(dataStates, ~ .x %in% availableStates)

  if (!is.null(dataStateUsed)) return(list(
    finalS4Object = getDIAStateReference(stateManager, dataStateUsed), dataStateUsed = dataStateUsed))
  catFn(sprintf("%s: No data S4 object found in any valid state\n", logPrefix))
  list(finalS4Object = NULL, dataStateUsed = NULL)
}

buildProtSummaryTemplateStatus <- function(projectDirs,
                                           omicType = "proteomics") {
  if (!omicType %in% names(projectDirs)) {
    return("[WARNING] Project directories not available")
  }

  templateDir <- file.path(
    projectDirs[[omicType]]$base_dir,
    "scripts",
    omicType
  )

  templateMap <- c(
    "DIA-NN" = "DIANN_report.rmd",
    "DIA-NN limpa" = "DIANN_limpa_report.rmd",
    "TMT" = "TMT_report.rmd",
    "LFQ" = "LFQ_report.rmd"
  )
  templatesStatus <- vapply(names(templateMap), function(label) {
    templatePath <- file.path(templateDir, templateMap[[label]])
    if (file.exists(templatePath)) {
      paste(label, "[OK]")
    } else {
      paste(label, "[missing]")
    }
  }, character(1))

  if (any(grepl("\\[OK\\]", templatesStatus))) {
    paste("Templates:", paste(templatesStatus, collapse = ", "))
  } else {
    "[WARNING] Report templates will be downloaded when generating report"
  }
}

resolveProtSummaryExpectedTemplate <- function(workflowType,
                                               limpaRequested = FALSE) {
  if (is.null(workflowType) || !nzchar(as.character(workflowType))) {
    workflowType <- "DIA"
  }
  workflowType <- tolower(as.character(workflowType))

  # KIT-046: the limpa report belongs to the DIA_LIMPA workflow type (`limpaRequested` is kept for callers that
  # already resolved the type to that flag).
  if (workflowType == "dia_limpa" || isTRUE(limpaRequested)) {
    return("DIANN_limpa_report.rmd")
  }

  if (workflowType %in% c("tmt", "tmt_pd")) {
    return("TMT_report.rmd")
  }

  if (workflowType == "lfq") {
    return("LFQ_report.rmd")
  }

  "DIANN_report.rmd"
}

resolveProtSummaryTemplateSelection <- function(workflowType,
                                                requestedTemplate = NULL,
                                                limpaRequested = FALSE) {
  expectedTemplate <- resolveProtSummaryExpectedTemplate(
    workflowType = workflowType,
    limpaRequested = limpaRequested
  )
  requestedTemplate <- if (!is.null(requestedTemplate) && nzchar(requestedTemplate)) {
    basename(requestedTemplate)
  } else {
    NULL
  }

  if (is.null(workflowType) || !nzchar(as.character(workflowType))) {
    workflowType <- "DIA"
  }
  compatibleTemplates <- switch(tolower(as.character(workflowType)),
    tmt = "TMT_report.rmd",
    tmt_pd = "TMT_report.rmd",
    lfq = "LFQ_report.rmd",
    dia = if (isTRUE(limpaRequested)) "DIANN_limpa_report.rmd" else "DIANN_report.rmd",
    dia_limpa = "DIANN_limpa_report.rmd",
    expectedTemplate
  )

  if (!is.null(requestedTemplate) && requestedTemplate %in% compatibleTemplates) {
    return(list(
      templateFilename = requestedTemplate,
      expectedTemplate = expectedTemplate,
      requestedTemplate = requestedTemplate,
      requestHonoured = TRUE,
      staleTemplateIgnored = FALSE
    ))
  }

  list(
    templateFilename = expectedTemplate,
    expectedTemplate = expectedTemplate,
    requestedTemplate = requestedTemplate,
    requestHonoured = FALSE,
    staleTemplateIgnored = !is.null(requestedTemplate)
  )
}

resolveProtSummaryReportTemplate <- function(workflowData,
                                             logPrefix = "REPORT",
                                             catFn = cat) {
  workflowTypeDetected <- NULL
  dataStateUsed <- NULL
  reportTemplateRequested <- NULL
  limpaRequested <- FALSE

  if (!is.null(workflowData) &&
      !is.null(workflowData$config_list) &&
      !is.null(workflowData$config_list$globalParameters)) {
    workflowTypeDetected <- workflowData$config_list$globalParameters$workflow_type
    reportTemplateRequested <- workflowData$config_list$globalParameters$report_template
    if (!is.null(reportTemplateRequested) && nzchar(reportTemplateRequested)) {
      reportTemplateRequested <- basename(reportTemplateRequested)
      catFn(sprintf(
        "   %s: Detected report_template from workflow_data: %s\n",
        logPrefix,
        reportTemplateRequested
      ))
    }
  }

  if (!is.null(workflowTypeDetected) && nzchar(workflowTypeDetected)) {
    catFn(sprintf(
      "   %s: Detected workflow_type from workflow_data: %s\n",
      logPrefix,
      workflowTypeDetected
    ))
  }

  if ((is.null(workflowTypeDetected) || !nzchar(workflowTypeDetected)) &&
      !is.null(workflowData) &&
      !is.null(workflowData$state_manager)) {
    finalS4State <- resolveProtSummaryFinalS4State(
      workflowData,
      logPrefix = logPrefix,
      catFn = catFn
    )
    currentS4 <- finalS4State$finalS4Object
    dataStateUsed <- finalS4State$dataStateUsed
    summaryArgs <- if (isSavedWorkflowReference(currentS4))
      readDIAStateMetadata(workflowData$state_manager, currentS4)$data$slots$args else
      tryCatch(currentS4@args, error = function(e) NULL)
    workflowTypeFromS4 <- tryCatch(
      summaryArgs$globalParameters$workflow_type,
      error = function(e) NULL
    )
    reportTemplateFromS4 <- tryCatch(
      summaryArgs$globalParameters$report_template,
      error = function(e) NULL
    )

    if ((is.null(reportTemplateRequested) || !nzchar(reportTemplateRequested)) &&
        !is.null(reportTemplateFromS4) &&
        nzchar(reportTemplateFromS4)) {
      reportTemplateRequested <- basename(reportTemplateFromS4)
      catFn(sprintf(
        "   %s: Detected report_template from S4 object (state: %s): %s\n",
        logPrefix,
        dataStateUsed,
        reportTemplateRequested
      ))
    }
    if (!is.null(workflowTypeFromS4) && nzchar(workflowTypeFromS4)) {
      workflowTypeDetected <- workflowTypeFromS4
      catFn(sprintf(
        "   %s: Detected workflow_type from S4 object (state: %s): %s\n",
        logPrefix,
        dataStateUsed,
        workflowTypeDetected
      ))
    }
  }

  # The workflow bus is the only configuration source (KIT-013); the global config_list is gone.
  if (is.null(workflowTypeDetected) || !nzchar(workflowTypeDetected)) {
    workflowTypeDetected <- "DIA"
    catFn(sprintf(
      "   %s: Using default workflow_type: DIA (no workflow type found)\n",
      logPrefix
    ))
  }
  # KIT-046: the limpa report follows the workflow type alone (never a flag or the presence of limpa arguments).
  limpaRequested <- identical(toupper(as.character(workflowTypeDetected)), "DIA_LIMPA")

  templateSelection <- resolveProtSummaryTemplateSelection(
    workflowType = workflowTypeDetected,
    requestedTemplate = reportTemplateRequested,
    limpaRequested = limpaRequested
  )
  templateFilename <- templateSelection$templateFilename

  if (isTRUE(templateSelection$staleTemplateIgnored)) {
    catFn(sprintf(
      "   %s: Ignoring stale report_template '%s' for workflow_type: %s; using %s\n",
      logPrefix,
      templateSelection$requestedTemplate,
      workflowTypeDetected,
      templateFilename
    ))
  }

  catFn(sprintf(
    "   %s: Selected template: %s for workflow_type: %s\n",
    logPrefix,
    templateFilename,
    workflowTypeDetected
  ))

  list(
    workflowTypeDetected = workflowTypeDetected,
    templateFilename = templateFilename,
    dataStateUsed = dataStateUsed,
    requestedTemplate = templateSelection$requestedTemplate,
    staleTemplateIgnored = templateSelection$staleTemplateIgnored
  )
}

summariseProtSummaryWorkflowMetadata <- function(workflowData,
                                                resolveReportTemplateFn = resolveProtSummaryReportTemplate) {
  globalParameters <- tryCatch(
    workflowData$config_list$globalParameters,
    error = function(e) NULL
  )
  workflowType <- globalParameters$workflow_type
  reportTemplate <- globalParameters$report_template

  if (!is.null(workflowData)) {
    templateInfo <- tryCatch(
      resolveReportTemplateFn(
        workflowData = workflowData,
        catFn = function(...) invisible(NULL)
      ),
      error = function(e) NULL
    )
    if (!is.null(templateInfo)) {
      workflowType <- templateInfo$workflowTypeDetected
      reportTemplate <- templateInfo$templateFilename
    }
  }

  list(
    workflow_type = workflowType,
    report_template = reportTemplate,
    parameters = list(
      globalParameters = globalParameters,
      da_ui_params = tryCatch(workflowData$da_ui_params, error = function(e) NULL),
      enrichment_ui_params = tryCatch(workflowData$enrichment_ui_params, error = function(e) NULL),
      # KIT-018: the accepted normalisation product first; the module copy for in-memory owners.
      ruv_optimization_result = tryCatch(resolveDIARuvOptimizationResult(workflowData$state_manager), error = function(e) NULL) %||%
        tryCatch(workflowData$ruv_optimization_result, error = function(e) NULL)
    )
  )
}

prepareProtSummarySessionStateExport <- function(projectDirs,
                                                 omicType = "proteomics",
                                                 experimentLabel,
                                                 description,
                                                 workflowArgsSaved,
                                                 filesCopied,
                                                 reportGenerated,
                                                 reportPath,
                                                 workflowData = NULL,
                                                 exportDate = Sys.Date(),
                                                 timestamp = Sys.time(),
                                                 summariseWorkflowMetadataFn = summariseProtSummaryWorkflowMetadata) {
  sessionExportPath <- file.path(
    projectDirs[[omicType]]$source_dir,
    paste0("session_state_", exportDate, ".RDS")
  )
  workflowMetadata <- summariseWorkflowMetadataFn(workflowData)

  sessionState <- list(
    experiment_label = experimentLabel,
    description = description,
    timestamp = timestamp,
    omic_type = omicType,
    workflow_args_saved = workflowArgsSaved,
    files_copied = filesCopied,
    report_generated = reportGenerated,
    report_path = reportPath,
    workflow_type = workflowMetadata$workflow_type,
    report_template = workflowMetadata$report_template,
    parameter_payload = workflowMetadata$parameters,
    project_dirs = projectDirs
  )

  list(sessionExportPath = sessionExportPath, sessionState = sessionState)
}
