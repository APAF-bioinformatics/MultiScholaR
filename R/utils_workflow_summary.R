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

# KIT-108 W3 (D3, D4): the Summary & Export tab's core (kit-101 design 4.10). Each entry takes an omic's `summary` group
# (`omicGuiConfig(<omic>)$summary`) and returns plain values; it reports through `notify(message, type, duration)` and
# `progress(value, detail)` callbacks and calls no Shiny UI function (the frame, mod_omic_summary.R, renders). The flow
# is the lipid, metab and trx handlers' at 5be34ee5 (mod_lipid_summary_{workflow,state}_helpers.R,
# mod_metab_summary_observer_helpers.R, mod_trx_summary.R; deleted by W3), with their per-omic differences read from the
# group: the final state, the record's owner class, the copy's required artifacts and `project_dirs` lookup, the
# integration file, its managed operation, its in-process writer and its gate, the report's template, directory and
# flow order, the session-state payload and the texts. The functions a test mocks (`copyToResultsSummary`,
# `RenderReport`, `createWorkflowArgsFromConfig`, `omicsManagedProduct`, `writeWorkflowAnalysisRecord`,
# `exportWorkflowAnalysisRecordBundle`, `utils::download.file`) are looked up by name when called. The copy and the
# render stay in process (D3). Transcriptomics' RUVIII export is an export step run by `runWorkflowSummaryExportStep`
# (D4). F108-03 (the lead's decision): a save whose integration write is refused by the gate shows no "Saved
# Integration" notice. The lipid, metab and trx functions below the entries are moved unchanged (D3's map).

# ------------------------------------------------------------------------------------------------
# The entries.
# ------------------------------------------------------------------------------------------------

#' A per-omic text of the group (`texts`), or the shared default.
#' @noRd
workflowSummaryText <- function(config, name, default) {
    value <- unname(config$texts[name])
    if (length(value) != 1L || is.na(value)) default else value
}

#' The summary's final state (D1's `final_state`): the omic's `context` (a reference, no object), or the first of its
#' data states found, head first, read as ids on a managed session (`read_ids`), else the state itself.
#' @param omicType The omic, for its managed runner.
#' @return `list(finalS4Object, contrastsTbl, dataStateUsed)`.
#' @noRd
workflowSummaryFinalState <- function(config, workflowData, omicType, catFn = cat) {
    final <- config$final_state
    if (!is.null(final$context)) {
        state <- final$context(workflowData, catFn = catFn)
        # KIT-108 W4: proteomics' context (`resolveProtSummaryFinalS4State`) names no contrasts; its save reads the bus's.
        if (!"contrastsTbl" %in% names(state)) state["contrastsTbl"] <- list(workflowData$contrasts_tbl)
        return(state)
    }
    finalS4Object <- NULL
    dataStateUsed <- NULL
    manager <- workflowData$state_manager
    if (!is.null(manager)) {
        # KIT-071: the final state is the head's (the owner's current state when it is a data state).
        dataStateUsed <- workflowSummaryFinalStateName(manager, final$states, detectFn = purrr::detect)
        if (!is.null(dataStateUsed)) {
            # KIT-082: on a managed session the worker writes the integration object, and the parameters read the ids.
            finalS4Object <- if (!is.null(workflowOmicsManagedActions(omicType))) final$read_ids(manager, dataStateUsed) else
                manager$getState(dataStateUsed)
            catFn(sprintf("SESSION SUMMARY: Retrieved DATA S4 object from state '%s'\n", dataStateUsed))
        } else {
            catFn("SESSION SUMMARY: No data S4 object found in any valid state\n")
        }
    } else {
        catFn("SESSION SUMMARY: No state manager available\n")
    }
    list(finalS4Object = finalS4Object, contrastsTbl = workflowData$contrasts_tbl, dataStateUsed = dataStateUsed)
}

#' Whether the analysis record (and its bundle) is written for this bus: always with a bus, or only for the group's
#' `record_class` owner.
#' @noRd
workflowSummaryRecords <- function(config, workflowData) {
    !is.null(workflowData) && (is.null(config$record_class) || inherits(workflowData$state_manager, config$record_class))
}

#' Save the workflow's parameters: the final state, `study_parameters.txt` (`createWorkflowArgsFromConfig`), the analysis
#' record, then the integration export. On an error the fallback `study_parameters.txt` is written (when absent).
#' @param label,description The experiment label and description.
#' @param notify `function(message, type, duration)`.
#' @param done `function(saved)` after the integration export: `saved` is `list(studyParamsFile, integrationFile)`.
#' @param failed `function(error)` after the fallback.
#' @noRd
saveWorkflowSummaryParameters <- function(config, workflowData, projectDirs, omicType, label, description, notify, done,
    failed, createWorkflowArgsFn = createWorkflowArgsFromConfig, writeAnalysisRecordFn = writeWorkflowAnalysisRecord,
    saveRDSFn = saveRDS, catFn = cat) {
    paths <- projectDirs[[omicType]]
    catFn("SESSION SUMMARY: Starting workflow args save process\n")
    summaryFailed <- function(error) {
        catFn("SESSION SUMMARY ERROR:", conditionMessage(error), "\n")
        basicParamsFile <- file.path(paths$source_dir, "study_parameters.txt")
        if (!file.exists(basicParamsFile)) {
            writeLines(c("Study Parameters", "================", "", paste("Workflow Name:", label),
                paste("Description:", description), paste("Timestamp:", Sys.time()),
                paste("Error:", conditionMessage(error))), basicParamsFile)
        }
        notify("Study parameters saved with warnings", "warning", 5)
        failed(error)
        invisible(basicParamsFile)
    }
    tryCatch({
        state <- workflowSummaryFinalState(config, workflowData, omicType, catFn)
        studyParamsFile <- createWorkflowArgsFn(workflow_name = label, description = description,
            source_dir_path = paths$source_dir, final_s4_object = state$finalS4Object, contrasts_tbl = state$contrastsTbl,
            workflow_data = workflowData)
        # KIT-031: the analysis record, derived from the owner's head, beside study_parameters.txt.
        if (workflowSummaryRecords(config, workflowData)) {
            recordFile <- file.path(paths$source_dir, "analysis_record.json")
            tryCatch(writeAnalysisRecordFn(workflowData$state_manager, recordFile), error = function(e) {
                notify(paste("The analysis record (analysis_record.json) was not written:", conditionMessage(e)), "warning",
                    NULL)
            })
        }
        integrationDir <- paths$integration_dir %||% file.path(paths$base_dir, "integration")
        if (!dir.exists(integrationDir)) dir.create(integrationDir, recursive = TRUE, showWarnings = FALSE)
        filename <- config$integration$filename
        filename <- if (is.function(filename)) filename(workflowData, omicType, label) else
            if (grepl("%", filename, fixed = TRUE)) sprintf(filename, omicType, label) else filename
        object <- state$finalS4Object
        stateName <- state$dataStateUsed
        state <- NULL # KIT-075: the object leaves this frame once written
        startWorkflowSummaryIntegrationExport(config, workflowData, omicType, file.path(integrationDir, filename),
            stateName, object, notify = notify, saveRDSFn = saveRDSFn, failed = summaryFailed,
            sourceDir = paths$source_dir, done = function(written) {
                # F108-03: a write the gate refused shows no "Saved Integration" notice. KIT-108 W4: a route writer's
                # notice is its own (proteomics': the native archive's, or the S4 RDS's).
                if (!isFALSE(written)) {
                    notify(if (is.character(written)) written else
                        workflowSummaryText(config, "integration_saved", "Saved Integration S4 Object"), "message", 5)
                }
                notify("Study parameters saved successfully", "message", 5)
                done(list(studyParamsFile = studyParamsFile, integrationFile = if (!isFALSE(written)) filename else
                    workflowSummaryText(config, "integration_not_written", filename)))
            })
    }, error = summaryFailed)
}

#' The integration export: refused (with the group's message) while its required export is not accepted on the head;
#' else the omic's managed product on its runner (the worker's `computeWorkflowOmicsIntegrationExportAction`), else the
#' in-process writer (`integration$write`, or the final object's RDS).
#' @param stateName The final state's name (the managed product's).
#' @param object The final object (the in-process RDS's content).
#' @param done `function(written)`: `written` is FALSE when the gate refused the write, else TRUE, or a route writer's
#'   notice.
#' @param sourceDir The project's source directory, for a route writer (KIT-108 W4: proteomics' native archive).
#' @noRd
startWorkflowSummaryIntegrationExport <- function(config, workflowData, omicType, path, stateName, object, notify, done,
    failed, saveRDSFn = saveRDS, sourceDir = NULL) {
    integration <- config$integration
    requires <- integration$requires
    if (!is.null(requires) && is.null(workflowData$state_manager$.context$snapshot$exports[[requires]])) {
        notify(integration$missing_message, "warning", 10)
        return(done(FALSE))
    }
    # KIT-082: the integration object is the whole state, so a managed session's worker writes it.
    runner <- if (!is.null(integration$operation)) workflowOmicsManagedActions(omicType)
    if (!is.null(runner) && !is.null(stateName)) {
        runner$run(list(omicsManagedProduct(integration$operation, list(path = path, state_name = stateName),
            function(result) invisible(NULL))), done = function() done(TRUE), failed = failed)
        return(invisible(path))
    }
    write <- integration$write
    # KIT-108 W4: a route writer (proteomics' `protSummaryIntegrationWrite`) also takes the final object, the source
    # directory and the RDS writer, and returns the notice the save shows.
    route <- !is.null(write) && "sourceDir" %in% names(formals(write))
    notice <- if (route) write(workflowData, path, object = object, sourceDir = sourceDir, saveRDSFn = saveRDSFn) else
        if (is.null(write)) saveRDSFn(object, path) else write(workflowData, path)
    object <- NULL
    done(if (route) notice else TRUE)
}

#' Copy the results to the publication directory: the fallback `study_parameters.txt` when the parameters were not
#' saved; the contrasts and design from the bus, else the source directory's tables; the `project_dirs` global, assigned
#' only when absent (`copy$inherits`); the gate (`copy$gate`); the exported `copyToResultsSummary(..., force = TRUE)`;
#' its failures classified by `copy$required` (a required one stops the copy).
#' @param argsSaved Whether the parameters were saved.
#' @return `list(fallbackWritten, required, optional, contrastsTbl, designMatrix, error)`; `error` is the failure's
#'   message, or NULL.
#' @noRd
copyWorkflowSummaryResults <- function(config, workflowData, projectDirs, omicType, label, description, argsSaved, notify,
    copyFn = copyToResultsSummary, catFn = cat) {
    paths <- projectDirs[[omicType]]
    fallbackWritten <- FALSE
    basicParamsFile <- file.path(paths$source_dir, "study_parameters.txt")
    if (!isTRUE(argsSaved) && !file.exists(basicParamsFile)) {
        fallbackWritten <- tryCatch({
            writeLines(paste("Study Parameters", "================", "", paste("Workflow Name:", label),
                paste("Description:", description), paste("Timestamp:", Sys.time()),
                "Note: Some parameters could not be saved due to serialization issues", sep = "\n"), basicParamsFile)
            TRUE
        }, error = function(e) {
            logger::log_error(logger::skip_formatter(paste("Failed to create basic study_parameters.txt:", conditionMessage(e))))
            FALSE
        })
    }
    tryCatch({
        contrastsTbl <- workflowData$contrasts_tbl
        designMatrix <- workflowData$design_matrix
        if (is.null(designMatrix) && file.exists(file.path(paths$source_dir, "design_matrix.tab"))) {
            designMatrix <- readr::read_tsv(file.path(paths$source_dir, "design_matrix.tab"), show_col_types = FALSE)
        }
        if (is.null(contrastsTbl) && file.exists(file.path(paths$source_dir, "contrasts_tbl.tab"))) {
            contrastsTbl <- readr::read_tsv(file.path(paths$source_dir, "contrasts_tbl.tab"), show_col_types = FALSE)
        }
        if (!is.null(config$copy$gate)) config$copy$gate()
        if (!exists("project_dirs", envir = .GlobalEnv, inherits = config$copy$inherits)) {
            catFn("SESSION SUMMARY: Setting project_dirs in global environment\n")
            assign("project_dirs", projectDirs, envir = .GlobalEnv)
        }
        failures <- copyFn(omic_type = omicType, experiment_label = label, contrasts_tbl = contrastsTbl,
            design_matrix = designMatrix, force = TRUE)
        required <- list()
        optional <- list()
        if (!is.null(config$copy$required) && length(failures) > 0L && !isTRUE(failures)) {
            isRequired <- vapply(failures, function(failure) {
                !is.null(failure$display_name) && failure$display_name %in% config$copy$required
            }, logical(1))
            required <- failures[isRequired]
            optional <- failures[!isRequired]
        }
        if (length(required) > 0L) {
            stop(sprintf("Required publication artifacts missing: %s",
                paste(vapply(required, function(failure) failure$display_name, character(1)), collapse = ", ")))
        }
        notify(if (length(optional) > 0L) "Publication files copied with optional missing artifacts" else
            "Publication files copied", if (length(optional) > 0L) "warning" else "message", 5)
        list(fallbackWritten = fallbackWritten, required = required, optional = optional, contrastsTbl = contrastsTbl,
            designMatrix = designMatrix, error = NULL)
    }, error = function(e) {
        notify(paste("Copy error:", conditionMessage(e)), "error", 10)
        logger::log_error(logger::skip_formatter(paste("Failed to copy files:", conditionMessage(e))))
        list(fallbackWritten = fallbackWritten, required = list(), optional = list(), contrastsTbl = NULL,
            designMatrix = NULL, error = conditionMessage(e))
    })
}

#' Generate the report: the template (`report$template`, a file name or the route's resolver), ensured in the project's
#' `scripts/<omic>/` (copied from the package's `reports/<dir>/`, else downloaded from GitHub), then the exported
#' `RenderReport`, then the analysis record bundle (before or after the success notice, `report$record_after_success`).
#' A failed template retrieval stops the report unless `report$render_after_failed_template`.
#' @param progress `function(value, detail)`.
#' @param within `function(message, code)`: runs the report under the frame's progress display.
#' @return `list(path, template, templateSource, bundleError)`; `path` is the rendered report, or NULL.
#' @noRd
renderWorkflowSummaryReport <- function(config, workflowData, projectDirs, omicType, label, notify, progress,
    within = function(message, code) code, renderFn = RenderReport, exportRecordFn = exportWorkflowAnalysisRecordBundle,
    downloadFileFn = utils::download.file, catFn = cat) {
    paths <- projectDirs[[omicType]]
    if (!omicType %in% names(projectDirs) || is.null(paths$base_dir)) {
        notify("Error: Project directories not properly initialized", "error", 10)
        return(invisible(list(path = NULL, template = NULL, templateSource = NULL, bundleError = NULL)))
    }
    report <- config$report
    within("Generating report...", {
        # KIT-108 W4: proteomics shows the route's detection first; its resolver returns the template with the route's
        # details (`templateFilename`).
        detect <- workflowSummaryText(config, "progress_detect", NA_character_)
        if (!is.na(detect)) progress(0.1, detect)
        template <- if (is.function(report$template)) report$template(workflowData) else report$template
        if (is.list(template)) template <- template$templateFilename
        catFn(sprintf("REPORT: Selected template: %s\n", template))
        templatePath <- file.path(paths$base_dir, "scripts", omicType, template)
        templateSource <- "existing"
        missing <- !file.exists(templatePath)
        if (missing || isTRUE(report$checking_progress_always)) {
            progress(0.2, paste("Checking for", template, "template..."))
        }
        if (missing) {
            dir.create(dirname(templatePath), recursive = TRUE, showWarnings = FALSE)
            templateSource <- tryCatch({
                packageFile <- system.file("reports", report$dir, template, package = "MultiScholaR")
                if (file.exists(packageFile) && packageFile != "") {
                    file.copy(packageFile, templatePath)
                    notify(paste(template, "template copied from package"), "message", 5)
                    "package"
                } else {
                    url <- paste0("https://raw.githubusercontent.com/APAF-bioinformatics/MultiScholaR/main/inst/reports/",
                        report$dir, "/", template)
                    downloadFileFn(url, destfile = templatePath, quiet = TRUE)
                    notify(paste(template, "template downloaded"), "message", 5)
                    "github"
                }
            }, error = function(e) {
                notify(paste("Template retrieval failed:", conditionMessage(e)), "error", 10)
                NULL
            })
        }
        if (missing && is.null(templateSource) && !isTRUE(report$render_after_failed_template)) {
            list(path = NULL, template = template, templateSource = NULL, bundleError = NULL)
        } else {
            progress(0.5, "Rendering report...")
            workflowSummaryRenderReport(config, workflowData, paths, omicType, label, template, templateSource, notify,
                progress, renderFn, exportRecordFn)
        }
    })
}

#' The render and the record bundle (`renderWorkflowSummaryReport`'s second half).
#' @noRd
workflowSummaryRenderReport <- function(config, workflowData, paths, omicType, label, template, templateSource, notify,
    progress, renderFn, exportRecordFn) {
    bundleError <- NULL
    exportBundle <- function() {
        if (!workflowSummaryRecords(config, workflowData)) return(invisible(NULL))
        # KIT-031: the analysis record bundle (JSON, Word, Excel) beside the report; a failure is shown and leaves the
        # report as generated (the KIT-029 pattern).
        progress(0.9, "Exporting the analysis record...")
        tryCatch(exportRecordFn(workflowData$state_manager, paths$study_report_dir %||% paths$results_summary_dir),
            error = function(e) {
                bundleError <<- conditionMessage(e)
                notify(paste("The analysis record bundle was not exported:", conditionMessage(e)), "warning", NULL)
            })
    }
    tryCatch({
        rendered <- renderFn(omic_type = omicType, experiment_label = label, rmd_filename = template)
        if (!is.null(rendered) && file.exists(rendered)) {
            if (!isTRUE(config$report$record_after_success)) exportBundle()
            notify("Report generated successfully!", "message", 5)
            if (isTRUE(config$report$record_after_success)) exportBundle()
        } else {
            notify("Report generation failed - no output file created", "error", 10)
            rendered <- NULL
        }
        list(path = rendered, template = template, templateSource = templateSource, bundleError = bundleError)
    }, error = function(e) {
        notify(paste("Report generation failed:", conditionMessage(e)), "error", 15)
        notify("Debug info: Check R console for detailed error trace", "warning", 10)
        list(path = NULL, template = template, templateSource = templateSource, bundleError = bundleError)
    })
}

#' Export the session state: `session_state_<date>.RDS` in the source directory, from the omic's payload
#' (`session_state`).
#' @param values The tab's reactive values, as a list.
#' @return The file written, or NULL after a failure (shown).
#' @noRd
writeWorkflowSummarySessionState <- function(config, workflowData, projectDirs, omicType, label, description, values,
    notify) {
    tryCatch({
        path <- file.path(projectDirs[[omicType]]$source_dir, paste0("session_state_", Sys.Date(), ".RDS"))
        state <- config$session_state(inputValues = list(experiment_label = label, description = description),
            projectDirs = projectDirs, omicType = omicType, values = values, timestamp = Sys.time(),
            workflowData = workflowData)
        saveRDS(state, path)
        notify(paste("Session state exported to:", path), "message", 5)
        path
    }, error = function(e) {
        notify(paste("Export failed:", conditionMessage(e)), "error", 10)
        logger::log_error(logger::skip_formatter(paste("Failed to export session state:", conditionMessage(e))))
        NULL
    })
}

#' The report template's status: the omic's (`report$status`), else the one template's.
#' @noRd
workflowSummaryTemplateStatus <- function(config, projectDirs, omicType) {
    if (is.function(config$report$status)) return(config$report$status(projectDirs, omicType))
    if (!omicType %in% names(projectDirs)) return("[WARNING] Project directories not available")
    template <- file.path(projectDirs[[omicType]]$base_dir, "scripts", omicType, config$report$template)
    if (file.exists(template)) sprintf("Template: %s Report [OK]", config$report$label) else
        "[WARNING] Report template will be downloaded when generating report"
}

#' Run an export step (D4; trx's RUVIII export, `runTrxRuviiiExport` at 5be34ee5 made generic over the step). On the
#' omic's runner: the step's product, accepted on the head by `accept`, then its publications as products. In process:
#' `run_reference`, then each publication (`runWorkflowOmicsPublicationInProcess`).
#' @param step One of the group's `export_steps`.
#' @param choice The step's choice (the input's value).
#' @param done `function(shown)`: `shown` is `list(text = status$text(value), at)`.
#' @param failed `function(error)`.
#' @param binding The omic's declaration binding.
#' @noRd
runWorkflowSummaryExportStep <- function(step, workflowData, experimentPaths, choice, done, failed, binding) {
    manager <- workflowData$state_manager
    runner <- workflowOmicsManagedActions(binding$omic_type)
    if (!is.null(runner)) {
        held <- new.env(parent = emptyenv())
        source <- getWorkflowStateReference(manager)
        export <- omicsManagedProduct(step$operation, step$parameters(manager, source, choice), function(result) {
            step$accept(manager, source, list(reference = result$reference, value = result$value))
            held$value <- result$value
        })
        actions <- c(list(export), lapply(step$publications, function(publication) {
            omicsManagedProduct(publication$operation, publication$parameters, function(result) invisible(NULL))
        }))
        runner$run(actions, done = function() done(list(text = step$status$text(held$value), at = Sys.time())),
            failed = failed)
        return(invisible(TRUE))
    }
    tryCatch({
        run <- step$run_reference(manager, choice)
        for (publication in step$publications) {
            runWorkflowOmicsPublicationInProcess(binding, manager, publication$operation, publication$parameters,
                experimentPaths)
        }
        done(list(text = step$status$text(run$value), at = Sys.time()))
    }, error = failed)
}

# ------------------------------------------------------------------------------------------------
# Transcriptomics' functions (moved unchanged from mod_trx_summary.R, D3; and the export step's parameters, D4).
# ------------------------------------------------------------------------------------------------

#' The M sources the export offers: the group, and the technical replicates only when the design has them.
#' @noRd
trxRuviiiMChoices <- function(has_replicates) {
    choices <- c(`Group (the notebook)` = "group")
    if (isTRUE(has_replicates)) choices <- c(choices, `Technical replicates` = "technical_replicate")
    choices
}

#' Whether the head's design has technical replicates (a technical-replicate set with more than one run); metadata only.
#' @noRd
trxHeadHasTechnicalReplicates <- function(manager) {
    meta <- tryCatch(readWorkflowOmicsMetadata(transcriptomicsDeclarationBinding, manager), error = function(error) NULL)
    column <- meta$slots$technical_replicate_id
    if (is.null(meta) || length(column) != 1L || is.na(column) || !column %in% names(meta$design_matrix)) return(FALSE)
    sets <- meta$design_matrix[[column]]
    any(duplicated(sets[!is.na(sets)]))
}

#' The status text of an accepted export.
#' @noRd
trxRuviiiStatusText <- function(value) {
    paste(c("RUVIII export complete (integration only; DA results unchanged).", sprintf("M = %s", value$m_source),
        sprintf("k = %d, %d inferred negative controls", as.integer(value$k), as.integer(value$control_count %||% sum(value$ctl_mask)))), collapse = "\n")
}

#' The summary's final state: the head's reference (its arguments are read from the record) and name; no object.
#' @noRd
trxSummaryContext <- function(workflowData, catFn = cat, ...) {
    manager <- workflowData$state_manager
    if (is.null(manager)) return(list(finalS4Object = NULL, contrastsTbl = workflowData$contrasts_tbl, dataStateUsed = NULL))
    catFn(sprintf("SESSION SUMMARY: final transcriptomics state '%s'\n", manager$current_state))
    list(finalS4Object = getWorkflowStateReference(manager), contrastsTbl = workflowData$contrasts_tbl, dataStateUsed = manager$current_state)
}

#' The head's metadata for the session export.
#' @noRd
trxSummaryStateObject <- function(workflowData) {
    tryCatch(readWorkflowOmicsMetadata(transcriptomicsDeclarationBinding, workflowData$state_manager), error = function(error) NULL)
}

#' The integration export in process: KIT-095's writer of the accepted RUVIII matrix.
#' @noRd
trxSummaryIntegrationWrite <- function(workflowData, path) {
    manager <- workflowData$state_manager
    withWorkflowStateScope(manager, function() trxWriteIntegrationExport(manager, list(), path))
}

#' Whether the RUVIII export the integration file is written from has been accepted.
#' @noRd
trxSummaryIntegrationReady <- function(workflowData) !is.null(workflowData$state_manager$.context$snapshot$exports$trx_ruviii_export)

#' The RUVIII export's parameters (D4; `runTrxRuviiiExport`'s at 5be34ee5, mod_trx_summary.R:76-78): the head's RUVs
#' logCPM and grid entities, and M.
#' @param source The head's reference.
#' @noRd
trxRuviiiExportParameters <- function(manager, source, m_source) {
    ruvs <- withWorkflowStateScope(manager, function() trxStateRuvs(manager, source))
    list(logcpm = ruvs$logcpm$entity_id, grid = ruvs$grid$entity_id, ruviii_m_source = m_source)
}

#' The session-state payload: lipid's, with the head's metadata (`trxSummaryStateObject`), the omic and its template
#' (`prepareLipidSummarySessionState` under the profile's entries at 5be34ee5, mod_transcriptomics.R:39, :45-57).
#' @noRd
trxSummarySessionState <- function(inputValues, projectDirs, omicType, values, timestamp, workflowData = NULL) {
    prepareLipidSummarySessionState(inputValues = inputValues, projectDirs = projectDirs, omicType = omicType,
        values = values, timestamp = timestamp, workflowData = workflowData, profile = list(omic_type = "transcriptomics",
            report_template = "transcriptomics_report.rmd", summary_state_object = trxSummaryStateObject))
}

# ------------------------------------------------------------------------------------------------
# Lipidomics' session-state functions (moved unchanged from mod_lipid_summary_state_helpers.R, D3).
# ------------------------------------------------------------------------------------------------

buildLipidSummarySessionState <- function(
    input,
    values,
    projectDirs,
    omicType,
    timeFn = Sys.time
) {
    list(
        experiment_label = input$experiment_label,
        description = input$description,
        timestamp = timeFn(),
        omic_type = omicType,
        workflow_args_saved = values$workflow_args_saved,
        files_copied = values$files_copied,
        report_generated = values$report_generated,
        report_path = values$report_path,
        project_dirs = projectDirs
    )
}

lipidSummaryFirstNonNull <- function(...) {
    values <- list(...)
    for (value in values) {
        if (!is.null(value)) {
            return(value)
        }
    }
    NULL
}

getLipidSummaryWorkflowValue <- function(workflowData, name) {
    if (is.null(workflowData)) {
        return(NULL)
    }

    tryCatch(workflowData[[name]], error = function(e) NULL)
}

getLipidSummaryStateObject <- function(workflowData) {
    stateManager <- getLipidSummaryWorkflowValue(workflowData, "state_manager")
    if (is.null(stateManager)) {
        return(NULL)
    }

    dataStates <- c(
        "lipid_correlation_filtered",
        "lipid_norm_complete",
        "lipid_ruv_corrected",
        "lipid_normalized",
        "loaded_for_de",
        "lipid_qc_complete"
    )

    # KIT-071: the final state is the head's (the owner's current state when it is a data state).
    stateName <- workflowSummaryFinalStateName(stateManager, dataStates)
    if (is.null(stateName) || is.na(stateName)) {
        stateName <- tryCatch(stateManager$current_state, error = function(e) NULL)
    }
    if (is.null(stateName) || is.na(stateName)) {
        return(NULL)
    }

    # KIT-075: the state's metadata (its slots, assay names and row counts, design), not the object.
    tryCatch(readLipidWorkflowMetadata(stateManager, stateName), error = function(e) NULL)
}

getLipidSummaryObjectArgs <- function(object) {
    if (is.list(object) && !isS4(object)) {
        return(object$slots$args) # KIT-075: the state's metadata
    }
    if (is.null(object) || !isS4(object) || !"args" %in% methods::slotNames(object)) {
        return(NULL)
    }

    tryCatch(object@args, error = function(e) NULL)
}

buildLipidSummaryParameterPayload <- function(workflowData, finalS4Object = NULL) {
    payload <- list(
        config_list = getLipidSummaryWorkflowValue(workflowData, "config_list"),
        contrasts_tbl = getLipidSummaryWorkflowValue(workflowData, "contrasts_tbl"),
        da_ui_params = getLipidSummaryWorkflowValue(workflowData, "da_ui_params"),
        normalization_ui_params = getLipidSummaryWorkflowValue(workflowData, "normalization_ui_params"),
        itsd_ui_params = getLipidSummaryWorkflowValue(workflowData, "itsd_ui_params"),
        ruv_optimization_result = getLipidSummaryWorkflowValue(workflowData, "ruv_optimization_result"),
        s4_args = getLipidSummaryObjectArgs(finalS4Object)
    )

    payload[!vapply(payload, is.null, logical(1))]
}

prepareLipidSummarySessionState <- function(
    inputValues,
    projectDirs,
    omicType,
    values,
    timestamp,
    workflowData = NULL,
    profile = NULL
) {
    sessionState <- list(
        experiment_label = inputValues$experiment_label,
        description = inputValues$description,
        timestamp = timestamp,
        omic_type = omicType,
        workflow_args_saved = values$workflow_args_saved,
        files_copied = values$files_copied,
        report_generated = values$report_generated,
        report_path = values$report_path,
        project_dirs = projectDirs
    )

    if (is.null(workflowData)) {
        return(sessionState)
    }

    finalS4Object <- if (is.null(profile)) getLipidSummaryStateObject(workflowData) else profile$summary_state_object(workflowData)
    parameterPayload <- buildLipidSummaryParameterPayload(workflowData, finalS4Object)
    globalParameters <- lipidSummaryFirstNonNull(
        parameterPayload$config_list$globalParameters,
        parameterPayload$s4_args$globalParameters,
        list()
    )

    sessionState$workflow_type <- lipidSummaryFirstNonNull(
        globalParameters$workflow_type,
        profile$omic_type %||% "lipidomics"
    )
    sessionState$report_template <- lipidSummaryFirstNonNull(
        globalParameters$report_template,
        profile$report_template %||% "lipidomics_report.rmd"
    )
    sessionState$parameter_payload <- parameterPayload

    if (is.list(finalS4Object) && !isS4(finalS4Object)) {
        # KIT-075: the state's metadata.
        sessionState$assay_names <- finalS4Object$assay_names
        sessionState$feature_counts <- finalS4Object$assay_rows
        sampleIdCol <- finalS4Object$slots$sample_id
        designMatrix <- finalS4Object$design_matrix
        if (is.data.frame(designMatrix) && length(sampleIdCol) == 1L && sampleIdCol %in% names(designMatrix)) {
            sessionState$sample_count <- length(unique(designMatrix[[sampleIdCol]]))
        }
    } else if (!is.null(finalS4Object) && isS4(finalS4Object)) {
        slotNames <- methods::slotNames(finalS4Object)
        if ("lipid_data" %in% slotNames) {
            assayData <- tryCatch(finalS4Object@lipid_data, error = function(e) NULL)
            if (is.list(assayData)) {
                sessionState$assay_names <- names(assayData)
                sessionState$feature_counts <- vapply(assayData, nrow, integer(1))
            }
        }
        if (all(c("design_matrix", "sample_id") %in% slotNames)) {
            designMatrix <- tryCatch(finalS4Object@design_matrix, error = function(e) NULL)
            sampleIdCol <- tryCatch(finalS4Object@sample_id, error = function(e) NULL)
            if (is.data.frame(designMatrix) && length(sampleIdCol) == 1L && sampleIdCol %in% names(designMatrix)) {
                sessionState$sample_count <- length(unique(designMatrix[[sampleIdCol]]))
            }
        }
    }

    contrastsTbl <- getLipidSummaryWorkflowValue(workflowData, "contrasts_tbl")
    if (is.data.frame(contrastsTbl)) {
        sessionState$contrast_count <- nrow(contrastsTbl)
    }

    sessionState
}

# ------------------------------------------------------------------------------------------------
# Metabolomics' session-state functions (moved unchanged from mod_metab_summary_observer_helpers.R, D3).
# ------------------------------------------------------------------------------------------------

metabSummaryFirstNonNull <- function(...) {
    values <- list(...)
    for (value in values) {
        if (!is.null(value)) {
            return(value)
        }
    }
    NULL
}

getMetabSummaryWorkflowValue <- function(workflowData, name) {
    if (is.null(workflowData)) {
        return(NULL)
    }

    tryCatch(workflowData[[name]], error = function(e) NULL)
}

getMetabSummaryStateObject <- function(workflowData) {
    state_manager <- getMetabSummaryWorkflowValue(workflowData, "state_manager")
    if (is.null(state_manager)) {
        return(NULL)
    }

    data_states <- c(
        "metab_correlation_filtered", "metab_norm_complete",
        "metab_ruv_corrected", "metab_normalized", "loaded_for_de",
        "metab_qc_complete"
    )

    # KIT-071: the final state is the head's (the owner's current state when it is a data state).
    state_name <- workflowSummaryFinalStateName(state_manager, data_states)
    if (is.null(state_name) || is.na(state_name)) {
        state_name <- tryCatch(state_manager$current_state, error = function(e) NULL)
    }
    if (is.null(state_name) || is.na(state_name)) {
        return(NULL)
    }

    # KIT-075: the state's metadata (its slots, assay names and row counts, design), not the object.
    tryCatch(readMetabWorkflowMetadata(state_manager, state_name), error = function(e) NULL)
}

getMetabSummaryObjectArgs <- function(object) {
    if (is.list(object) && !isS4(object)) {
        return(object$slots$args) # KIT-075: the state's metadata
    }
    if (is.null(object) || !isS4(object) || !"args" %in% methods::slotNames(object)) {
        return(NULL)
    }

    tryCatch(object@args, error = function(e) NULL)
}

buildMetabSummaryParameterPayload <- function(workflowData, finalS4Object = NULL) {
    payload <- list(
        config_list = getMetabSummaryWorkflowValue(workflowData, "config_list"),
        contrasts_tbl = getMetabSummaryWorkflowValue(workflowData, "contrasts_tbl"),
        da_ui_params = getMetabSummaryWorkflowValue(workflowData, "da_ui_params"),
        normalization_ui_params = getMetabSummaryWorkflowValue(workflowData, "normalization_ui_params"),
        itsd_ui_params = getMetabSummaryWorkflowValue(workflowData, "itsd_ui_params"),
        ruv_optimization_result = getMetabSummaryWorkflowValue(workflowData, "ruv_optimization_result"),
        enrichment_ui_params = getMetabSummaryWorkflowValue(workflowData, "enrichment_ui_params"),
        s4_args = getMetabSummaryObjectArgs(finalS4Object)
    )

    payload[!vapply(payload, is.null, logical(1))]
}

prepareMetabSummarySessionState <- function(
    inputValues,
    projectDirs,
    omicType,
    values,
    timestamp,
    workflowData = NULL
) {
    session_state <- list(
        experiment_label = inputValues$experiment_label,
        description = inputValues$description,
        timestamp = timestamp,
        omic_type = omicType,
        workflow_args_saved = values$workflow_args_saved,
        files_copied = values$files_copied,
        report_generated = values$report_generated,
        report_path = values$report_path,
        project_dirs = projectDirs
    )

    if (is.null(workflowData)) {
        return(session_state)
    }

    final_s4_object <- getMetabSummaryStateObject(workflowData)
    parameter_payload <- buildMetabSummaryParameterPayload(workflowData, final_s4_object)
    global_parameters <- metabSummaryFirstNonNull(
        parameter_payload$config_list$globalParameters,
        parameter_payload$s4_args$globalParameters,
        list()
    )

    session_state$workflow_type <- metabSummaryFirstNonNull(
        global_parameters$workflow_type,
        "metabolomics"
    )
    session_state$report_template <- metabSummaryFirstNonNull(
        global_parameters$report_template,
        "metabolomics_report.rmd"
    )
    session_state$parameter_payload <- parameter_payload

    if (is.list(final_s4_object) && !isS4(final_s4_object)) {
        # KIT-075: the state's metadata.
        session_state$assay_names <- final_s4_object$assay_names
        session_state$feature_counts <- final_s4_object$assay_rows
        sample_id_col <- final_s4_object$slots$sample_id
        design_matrix <- final_s4_object$design_matrix
        if (is.data.frame(design_matrix) && length(sample_id_col) == 1L && sample_id_col %in% names(design_matrix)) {
            session_state$sample_count <- length(unique(design_matrix[[sample_id_col]]))
        }
    } else if (!is.null(final_s4_object) && isS4(final_s4_object)) {
        slot_names <- methods::slotNames(final_s4_object)
        if ("metabolite_data" %in% slot_names) {
            assay_data <- tryCatch(final_s4_object@metabolite_data, error = function(e) NULL)
            if (is.list(assay_data)) {
                session_state$assay_names <- names(assay_data)
                session_state$feature_counts <- vapply(assay_data, nrow, integer(1))
            }
        }
        if (all(c("design_matrix", "sample_id") %in% slot_names)) {
            design_matrix <- tryCatch(final_s4_object@design_matrix, error = function(e) NULL)
            sample_id_col <- tryCatch(final_s4_object@sample_id, error = function(e) NULL)
            if (is.data.frame(design_matrix) && length(sample_id_col) == 1L && sample_id_col %in% names(design_matrix)) {
                session_state$sample_count <- length(unique(design_matrix[[sample_id_col]]))
            }
        }
    }

    contrasts_tbl <- getMetabSummaryWorkflowValue(workflowData, "contrasts_tbl")
    if (is.data.frame(contrasts_tbl)) {
        session_state$contrast_count <- nrow(contrasts_tbl)
    }

    session_state
}
