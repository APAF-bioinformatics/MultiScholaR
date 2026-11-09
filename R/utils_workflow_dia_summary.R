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
# resolves, as the report step does.

#' The copy's gate: the session's owner has no export pending (today's `assertDIAExportReady(getDIAActiveManager())`).
#' @noRd
protSummaryCopyGate <- function() assertDIAExportReady(getDIAActiveManager())

#' The integration file's name: the native archive on a native session, else the final object's RDS.
#' @noRd
protSummaryIntegrationFilename <- function(workflowData, omicType, label) {
  native <- identical(omicType, "proteomics") && isDIANativeSession(workflowData$state_manager)
  if (native) sprintf("%s_%s_final.msdia.tar", omicType, label) else sprintf("%s_%s_final_s4.RDS", omicType, label)
}

#' The integration file, written in process: on a native session the summary's archive is refreshed and copied (no
#' assay or S4 object is read); else the final object (a state reference) is saved as an RDS.
#' @param object The final object (`final_state$context`'s).
#' @param sourceDir The project's source directory (the native archive's).
#' @return The notice the save shows.
#' @noRd
protSummaryIntegrationWrite <- function(workflowData, path, object = NULL, sourceDir = NULL, saveRDSFn = saveRDS) {
  # The archive's name is `protSummaryIntegrationFilename`'s choice on a native session.
  native <- endsWith(path, ".msdia.tar")
  if (native) {
    refreshDIANativeSummaryArchive(workflowData, sourceDir)
    copyDIANativeSessionArtifact(sourceDir, path)
  } else {
    saveRDSFn(object, path)
  }
  if (native) "Saved Integration artifact" else "Saved Integration S4 Object"
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
