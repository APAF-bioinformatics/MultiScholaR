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

#' KIT-207 (VAL-105 (2), VAL-148): whether the Summary may write from the head (its final object, study parameters,
#' analysis record, integration file, report bundle and session state). Lipidomics and metabolomics accept the heads the
#' DA's head rule accepts (KIT-194: `omicDaFallbackHeadNormalised`; a head adopted from a session file unless
#' `omicDaLoadedStateRefused` names it), without the DA's precondition on the normalisation status; transcriptomics
#' accepts a normalised state (`trxDaReadMetadata`). Every other omic, and a manager that is not a session owner, as
#' before. Recorded divergence: the refusals use plain words for non-specialists, not conventions/R.md's "state the
#' offending value" (as KIT-194).
#' @param workflowData The workflow bus.
#' @param omicType The omic.
#' @return TRUE, invisibly; otherwise a condition of class `omic_norm_correlation_input_stale` whose text names the case.
#' @noRd
assertWorkflowSummaryHeadFinal <- function(workflowData, omicType) {
    manager <- workflowData$state_manager
    # KIT-216 (I10, VAL-158 parity): the proteomics rule by name, for a reference-state owner only.
    if (identical(omicType, "proteomics")) return(assertProtSummaryHeadFinal(manager))
    if (!isTRUE(omicType %in% c("lipidomics", "metabolomics", "transcriptomics")) || !inherits(manager, "DIAStateAdapter")) {
        return(invisible(TRUE))
    }
    binding <- omicGuiConfig(omicType)$binding$binding
    texts <- omicNormGuardTexts()
    refusal <- withWorkflowStateScope(manager, function() {
        if (omicNormIsTrx(binding)) {
            if (!is.null(trxDaReadMetadata(workflowData))) return(NULL)
            return(if (is.null(getWorkflowStateReference(manager, "trx_ruvs"))) texts$summary_none else texts$summary_trx)
        }
        head <- manager$current_state
        if (is.character(head) && length(head) == 1L && omicNormAdopted(manager, head)) {
            return(if (omicDaLoadedStateRefused(binding, head)) sprintf(texts$summary_loaded, head))
        }
        if (omicDaFallbackHeadNormalised(workflowData, binding)) return(NULL)
        if (is.null(getWorkflowStateReference(manager, paste0(binding$prefix, "_normalized")))) texts$summary_none else
            texts$summary
    })
    if (is.null(refusal)) return(invisible(TRUE))
    rlang::abort(refusal, class = "omic_norm_correlation_input_stale")
}

#' KIT-216 (I10; P05, P11; VAL-158, VAL-179 SD-1): the proteomics Summary's head rule, KIT-207's in proteomics words.
#' The head's table (the head; a result head, which only an unconverted store holds, by KIT-212's
#' `workflowDaTableState`) must be one of the proteomics DA's states, by name; no comparison with an export (KIT-206's
#' guard stays in the archive refresh). Refused with `summary_none` when no normalisation state is in the store,
#' otherwise with the proteomics text. Other owners (in-memory, test fakes, another omic's owner) pass.
#' @param manager The bus's state manager.
#' @return TRUE, invisibly; otherwise a condition of class `omic_norm_correlation_input_stale`.
#' @noRd
assertProtSummaryHeadFinal <- function(manager) {
    if (!usesDIAReferenceStates(manager)) return(invisible(TRUE))
    # r1.4 (audit F04): fail-safe; only an owner of another omic's binding (lipidomics, metabolomics, transcriptomics)
    # is skipped, so every proteomics owner of any route stays guarded.
    if (isTRUE(toupper(manager$workflow_type) %in% toupper(omicsWorkflowTypes()))) return(invisible(TRUE))
    texts <- omicNormGuardTexts()
    refusal <- withWorkflowStateScope(manager, function() {
        head <- manager$current_state
        table <- tryCatch({
            reference <- workflowDaTableState(manager)$reference
            snapshot <- manager$.context$snapshot
            found <- names(snapshot$states)[vapply(snapshot$states, function(state) identical(state$entity_id,
                reference$entity_id), logical(1))]
            if (head %in% found) head else utils::tail(intersect(diaStateHistory(manager), found), 1L)
        }, error = function(e) head)
        if (isTRUE(table %in% protSummaryDaStates)) return(NULL)
        normalised <- any(vapply(c("normalized", "ruv_corrected", "ruv_covariates"), function(state)
            !is.null(getDIAStateReference(manager, state)), logical(1)))
        if (normalised) texts$summary_prot else texts$summary_prot_none
    })
    if (is.null(refusal)) return(invisible(TRUE))
    rlang::abort(refusal, class = "omic_norm_correlation_input_stale")
}

#' KIT-212 (VAL-172 (1), (5)): whether the Summary may present the store's differential analysis: refused when it was run
#' on another table than the current one (`workflowDaResultStatus`), for every omic; with `files` (Copy and Report) also
#' when the omic's DA folders hold DA files no DA of this session made (`workflowDaFilesUnmade`), up to three named.
#' @param workflowData The workflow bus.
#' @param omicType The omic.
#' @param paths The omic's project paths.
#' @param files Whether the DA files on disk are checked.
#' @return TRUE, invisibly; otherwise a condition of class `workflow_summary_da_stale`.
#' @noRd
assertWorkflowSummaryDaCurrent <- function(workflowData, omicType, paths, files = FALSE) {
    manager <- workflowData$state_manager
    texts <- omicNormGuardTexts()
    if (identical(workflowDaResultStatus(manager), "stale")) rlang::abort(texts$summary_da_stale, class = "workflow_summary_da_stale")
    if (isTRUE(files)) {
        unmade <- workflowDaFilesUnmade(manager, paths, omicType)
        if (length(unmade)) rlang::abort(workflowDaUnmadeText(unmade, texts), class = "workflow_summary_da_stale")
    }
    invisible(TRUE)
}

#' KIT-243 (part D): refused while the newest DA's publication has not completed (its marker,
#' `removeWorkflowDaOwnFigures`, is still in the publication folder), for every omic.
#' @param paths The omic's project paths.
#' @return TRUE, invisibly; otherwise a condition of class `workflow_da_publication_pending`.
#' @noRd
assertWorkflowDaPublicationComplete <- function(paths) {
    marker <- workflowDaPublicationMarker(paths)
    if (!is.null(marker) && file.exists(marker)) {
        rlang::abort(omicNormGuardTexts()$da_publication_pending, class = "workflow_da_publication_pending")
    }
    invisible(TRUE)
}

#' KIT-243: what part B compares of an analysis record: its steps as (state, result entity) pairs in order, and its
#' differential-analysis block, after one JSON round trip with the record writer's options (`writeWorkflowAnalysisRecord`,
#' which is not called), so a fresh record and one read from disk are compared in the same form.
#' @param record An analysis record, built or read with `jsonlite::read_json(simplifyVector = FALSE)`.
#' @return A list.
#' @noRd
workflowSummarySaveIdentity <- function(record) {
    identity <- list(steps = lapply(record$steps, function(step) list(state = step$state, result_entity_id = step$result_entity_id)),
        differential_abundance = record$analyses$differential_abundance)
    jsonlite::fromJSON(jsonlite::toJSON(identity, auto_unbox = TRUE, digits = I(17), null = "null", na = "null"),
        simplifyVector = FALSE)
}

#' KIT-243 (VAL-212 (2)): whether the last "Save Workflow Arguments" describes the current table and analysis: the
#' on-disk `analysis_record.json`'s steps and DA block against the record the owner's store gives now
#' (`buildWorkflowAnalysisRecord`, the builder Save uses; no assay hydrated). "none" for a manager that is not a
#' reference-state session owner, when no record was saved, or when one cannot be read or compared (logged; KIT-212
#' E-5's fail-open rule, the recorded divergence over conventions/R.md's specific conditions).
#' @param manager The bus's state manager.
#' @param sourceDir The omic's source directory (where Save writes the record).
#' @return `"none"`, `"current"` or `"changed"`.
#' @noRd
workflowSummarySaveStatus <- function(manager, sourceDir) {
    if (!usesWorkflowReferenceStates(manager) || is.null(sourceDir)) return("none")
    path <- file.path(sourceDir, "analysis_record.json")
    if (!file.exists(path)) return("none")
    tryCatch({
        saved <- workflowSummarySaveIdentity(jsonlite::read_json(path, simplifyVector = FALSE))
        current <- workflowSummarySaveIdentity(withWorkflowStateScope(manager, function() buildWorkflowAnalysisRecord(manager)))
        if (identical(saved, current)) "current" else "changed"
    }, error = function(error) {
        logger::log_warn(logger::skip_formatter(paste("The saved analysis record could not be compared with the session:",
            conditionMessage(error))))
        "none"
    })
}

#' KIT-243 (VAL-212 (2)): Copy and the Report refuse until "Save Workflow Arguments" is pressed again when the last Save
#' does not describe the current table and analysis (`workflowSummarySaveStatus`). Save itself is not guarded.
#' @param workflowData The workflow bus.
#' @param paths The omic's project paths.
#' @return TRUE, invisibly; otherwise a condition of class `workflow_summary_save_stale`.
#' @noRd
assertWorkflowSummarySaveCurrent <- function(workflowData, paths) {
    if (identical(workflowSummarySaveStatus(workflowData$state_manager, paths$source_dir), "changed")) {
        rlang::abort(omicNormGuardTexts()$summary_save_stale, class = "workflow_summary_save_stale")
    }
    invisible(TRUE)
}

#' KIT-219 (VAL-187; drafts for the owner's wording pass): the texts of the enrichment files' notice and refusal.
#' @return A list of `sprintf` formats.
#' @noRd
workflowEnrichmentTexts <- function() {
    list(
        not_copied = paste("%d enrichment file(s) in the pathway_enrichment folder were not copied, because the current",
            "enrichment did not make them (they come from an earlier enrichment or an earlier analysis). They are still",
            "on disk: %s."),
        not_copied_no_analysis = paste("%d enrichment file(s) in the pathway_enrichment folder were not copied: this",
            "application has no enrichment analysis for %s, so they were not made by this analysis. They are still on",
            "disk: %s."),
        report_not_current = paste("The report was not generated: the publication folder holds enrichment results that",
            "the current enrichment did not make. Press 'Copy to Publication Directory' again, then generate the report."))
}

#' KIT-219: the sentence of the Copy notice about the enrichment files left out, up to three names.
#' @param record The copy's `enrichment_not_current` record (`files`).
#' @param omicType The omic.
#' @return A string.
#' @noRd
workflowEnrichmentNotCopiedText <- function(record, omicType) {
    texts <- workflowEnrichmentTexts()
    files <- record$files
    listed <- paste(utils::head(files, 3L), collapse = ", ")
    if (length(files) > 3L) listed <- sprintf("%s and %d more", listed, length(files) - 3L)
    if (identical(omicType, "proteomics")) sprintf(texts$not_copied, length(files), listed) else
        sprintf(texts$not_copied_no_analysis, length(files), omicType, listed)
}

#' The names one enrichment result's files have in `pathway_dir`: per contrast and direction with a result, the table with
#' the contrast as written (`<contrast>_<direction>_enrichment_results.tsv`) and the interactive and static plot with the
#' contrast sanitised as the kernel sanitises it (`<contrast>_<direction>_enrichment_plot.html` and `.png`). A direction
#' with a result but no terms writes no file; the names are those a run can write, never fewer.
#' @param results An `EnrichmentResults`.
#' @return Character vector of file names.
#' @noRd
workflowEnrichmentResultFileNames <- function(results) {
    data <- methods::slot(results, "enrichment_data")
    as.character(unlist(lapply(names(data), function(contrast) {
        directions <- intersect(c("up", "down"), names(Filter(Negate(is.null), data[[contrast]])))
        safe <- gsub("[^A-Za-z0-9_.-]", "_", contrast)
        lapply(directions, function(direction) c(paste0(contrast, "_", direction, "_enrichment_results.tsv"),
            paste0(safe, "_", direction, "_enrichment_plot.", c("html", "png"))))
    })))
}

#' KIT-219 (VAL-187): the files of `pathway_dir` that the current enrichment made, for Copy and the Report. A
#' reference-state session owner of proteomics: the files of the store's current enrichment product (`diaEnrichmentProduct`)
#' by `workflowEnrichmentResultFileNames`; `character()` when there is no current product, or when the record cannot be
#' read (logged; the divergence from `conventions/R.md`'s specific conditions is KIT-212 E-5's fail-open form, here the
#' safe side: nothing is presented). Another omic has no enrichment analysis, so none of its files is current:
#' `character()`. A manager that is not a reference-state session owner (a test fake, the legacy staging state) is not
#' checked: `NULL`.
#' @param workflowData The workflow bus.
#' @param omicType The omic.
#' @return `NULL` (unchecked), or a character vector of file names.
#' @noRd
workflowEnrichmentCurrentFiles <- function(workflowData, omicType) {
    manager <- workflowData$state_manager
    if (!usesWorkflowReferenceStates(manager)) return(NULL)
    if (!identical(omicType, "proteomics")) return(character())
    tryCatch(withWorkflowStateScope(manager, function() {
        marker <- diaEnrichmentProduct(manager)
        if (is.null(marker)) return(character())
        workflowEnrichmentResultFileNames(readDIAStateObject(manager, marker$state))
    }), error = function(error) {
        logger::log_warn(logger::skip_formatter(paste("The current enrichment's files could not be read:",
            conditionMessage(error))))
        character()
    })
}

#' KIT-219 (VAL-187): the Report is refused when the publication folder's `Enrichment_Plots` holds a file (listed
#' recursively) that the current enrichment did not make, or whose bytes differ from the same-named file in `pathway_dir`
#' (a missing source counts as different); "Copy again" then clears the folder and copies the current files only. A
#' second form beside KIT-212's refusal of DA files, stated as such: Copy omits and says, the Report refuses on what Copy
#' left in the folder. Not checked for a manager that is not a reference-state session owner.
#' @param workflowData The workflow bus.
#' @param omicType The omic.
#' @param paths The omic's project paths (`results_summary_dir`, `pathway_dir`).
#' @return TRUE, invisibly; otherwise a condition of class `workflow_summary_enrichment_not_current`.
#' @noRd
assertWorkflowSummaryEnrichmentCurrent <- function(workflowData, omicType, paths) {
    current <- workflowEnrichmentCurrentFiles(workflowData, omicType)
    if (is.null(current)) return(invisible(TRUE))
    folder <- file.path(paths$results_summary_dir %||% "", "Publication_figures", "Enrichment_Plots")
    shown <- if (dir.exists(folder)) list.files(folder, recursive = TRUE) else character()
    if (!length(shown)) return(invisible(TRUE))
    copied <- unname(tools::md5sum(file.path(folder, shown)))
    source <- unname(tools::md5sum(file.path(paths$pathway_dir %||% "", shown)))
    same <- shown %in% current & !is.na(copied) & !is.na(source) & copied == source
    if (!all(same)) {
        rlang::abort(workflowEnrichmentTexts()$report_not_current, class = "workflow_summary_enrichment_not_current")
    }
    invisible(TRUE)
}

#' The files refusal's text: the folder that holds every listed file, and up to three of them relative to it.
#' @param files Paths.
#' @param texts `omicNormGuardTexts()`.
#' @noRd
workflowDaUnmadeText <- function(files, texts) {
    folder <- dirname(files[[1L]])
    while (!all(startsWith(files, paste0(folder, "/"))) && !identical(dirname(folder), folder)) folder <- dirname(folder)
    shown <- substring(files, nchar(folder) + 2L)
    listed <- paste(utils::head(shown, 3L), collapse = ", ")
    if (length(shown) > 3L) listed <- sprintf("%s and %d more", listed, length(shown) - 3L)
    sprintf(texts$summary_da_files_unmade, folder, listed)
}

#' Whether the analysis record (and its bundle) is written for this bus: always with a bus, or only for the group's
#' `record_class` owner.
#' @noRd
workflowSummaryRecords <- function(config, workflowData) {
    !is.null(workflowData) && (is.null(config$record_class) || inherits(workflowData$state_manager, config$record_class))
}

#' KIT-228 (VAL-201; DRAFT wording, VAL-41): what the user is told when the omic's earlier integration input was removed
#' and this Save wrote no new one. `kind`: "failed" (an error or refusal, whose message is above it), "failed_unnamed" (an
#' error whose reason the Summary did not show), "gate" (the transcriptomics RUVIII export is missing) or "no_table" (a
#' writer or worker that wrote no input).
#' @param omicType The omic.
#' @param kind One of the three.
#' @return A string.
#' @noRd
workflowSummaryInputRemovedText <- function(omicType, kind) {
    folder <- sprintf("integration/%s/", omicType)
    switch(kind,
        gate = sprintf(paste("This omic's earlier integration input (%s) was removed because this save did not write a new",
            "one: the RUVIII export is missing. Run the RUVIII export, then save again before integrating."), folder),
        no_table = sprintf(paste("This omic's earlier integration input (%s) was removed, and this save did not write a new",
            "one. Save this Summary again before integrating; if this message returns, finish the analysis first (the",
            "integration input is written from the final normalised table)."), folder),
        failed_unnamed = sprintf(paste("This omic's earlier integration input (%s) was removed because this save did not write a",
            "new one. The integration cannot use %s until this Summary is saved without errors."), folder, omicType),
        sprintf(paste("This omic's earlier integration input (%s) was removed because this save did not write a new one (the",
            "reason is in the message above). The integration cannot use %s until this Summary is saved without errors."),
            folder, omicType))
}

#' KIT-228 (VAL-393; DRAFT wording, VAL-41): said once when a Save that removed the omic's earlier integration input is
#' handed to a runner that is retrying a failed publication.
#' @param omicType The omic.
#' @return A string.
#' @noRd
workflowSummaryInputRetryText <- function(omicType) {
    sprintf(paste("This omic's earlier integration input (integration/%s/) was removed by this save. If this save only finishes",
        "an earlier save's publication, that does not write it again: save once more before integrating."), omicType)
}

#' KIT-228 (DRAFT wording, VAL-41): the earlier integration input could not be removed at the start of the Save, and this
#' Save wrote no new one.
#' @param folder The omic's input folder.
#' @return A string.
#' @noRd
workflowSummaryInputNotRemovedText <- function(folder) {
    sprintf(paste("The earlier integration input at %s could not be removed. The integration may still read it, although",
        "this save did not write it: delete that folder, or save again once it can be removed."), folder)
}

#' Save the workflow's parameters: first the omic's earlier integration input is removed (KIT-228: one removal for every
#' way the Save can end without a new input, so the integration cannot import an input this Save did not write; a removal
#' that fails is named and the Save goes on), then the final state, `study_parameters.txt` (`createWorkflowArgsFromConfig`),
#' the analysis record and the integration export. On an error the fallback `study_parameters.txt` is written (when absent).
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
    # KIT-228: the omic's input folder (`integration/<omic>/`, where every writer puts it), whether one existed before this
    # Save and whether its removal failed.
    earlierDir <- paths$integration_dir %||% if (!is.null(paths$base_dir)) file.path(paths$base_dir, "integration")
    inputFolder <- if (!is.null(earlierDir)) file.path(earlierDir, omicType)
    earlier <- FALSE
    removalFailed <- FALSE
    # Whether a failure's reason was shown to the user: the managed step's own notice carries it (a failure from the runner).
    reasonShown <- FALSE
    # KIT-228 (E-5; audit I06): said only when an earlier input existed and this Save ends without a new one, chosen from
    # the state at the END of the Save: "could not be removed" only when the start removal failed and the folder is still
    # there, else that it was removed. `kind`: "failed" (an error or refusal), "failed_unnamed", "gate", "no_table".
    notifyInputRemoved <- function(kind) {
        if (!earlier) return(invisible(NULL))
        if (removalFailed && dir.exists(inputFolder)) {
            notify(workflowSummaryInputNotRemovedText(inputFolder), "warning", NULL)
        } else if (!file.exists(file.path(inputFolder, INTEGRATION_INPUT_MANIFEST))) {
            notify(workflowSummaryInputRemovedText(omicType, kind), "warning", NULL)
        }
        invisible(NULL)
    }
    summaryFailed <- function(error) {
        catFn("SESSION SUMMARY ERROR:", conditionMessage(error), "\n")
        # KIT-207: the head is not the normalisation's result. Refused before anything was written: no fallback file, no
        # "saved with warnings"; the omic's earlier integration input was removed at the start of the Save (SD-5 (b), K061-3;
        # KIT-228) and the user is told.
        # KIT-212 (VAL-172 (1)): the DA shown was run on another table; refused the same way.
        if (inherits(error, c("omic_norm_correlation_input_stale", "workflow_summary_da_stale"))) {
            notify(paste("The session could not be saved:", cli::ansi_strip(conditionMessage(error))), "error", NULL)
            notifyInputRemoved("failed")
            failed(error)
            return(invisible(NULL))
        }
        # KIT-052 (KIT-134 run 44379470655): a file the save could not write (the session archive's, on Windows with long
        # paths disabled) is named to the user, not only logged; the copy and the report stay refused until it is written.
        # KIT-052, KIT-206; KIT-207 (SD-4): a file the save could not write, or refused, is named to the user.
        named <- inherits(error, c("workflow_figure_unwritten", "prot_norm_correlation_input_stale", "workflow_export_dependency_missing",
            "integration_input_invalid"))
        if (named) {
            notify(paste("The session could not be saved:", cli::ansi_strip(conditionMessage(error))), "error", NULL)
        }
        basicParamsFile <- file.path(paths$source_dir, "study_parameters.txt")
        if (!file.exists(basicParamsFile)) {
            writeLines(c("Study Parameters", "================", "", paste("Workflow Name:", label),
                paste("Description:", description), paste("Timestamp:", Sys.time()),
                paste("Error:", conditionMessage(error))), basicParamsFile)
        }
        notify("Study parameters saved with warnings", "warning", 5)
        notifyInputRemoved(if (named || reasonShown) "failed" else "failed_unnamed")
        failed(error)
        invisible(basicParamsFile)
    }
    tryCatch({
        # KIT-228 (VAL-201): the omic's earlier integration input is removed once, here, before anything is written: every
        # Save then either publishes a new input or leaves none, whatever ends it (a refusal, an error, the gate, a worker
        # failure whose class is lost, an interrupted Save). A removal that fails is named and the Save goes on.
        # Also a Save pressed again while the runner retries a failed publication (VAL-393, audit C01 and C02): the runner
        # reports a retry for ANY step of the omic with a publication pending, and a resumed request can be dropped as
        # superseded, so the input is removed here as on every Save; the retry path says so (below).
        if (!is.null(inputFolder)) {
            earlier <- dir.exists(inputFolder)
            tryCatch(removeOmicIntegrationInput(earlierDir, omicType),
                error = function(e) notify(cli::ansi_strip(conditionMessage(e)), "warning", NULL))
            removalFailed <- dir.exists(inputFolder)
        }
        assertWorkflowSummaryHeadFinal(workflowData, omicType)
        assertWorkflowSummaryDaCurrent(workflowData, omicType, paths)
        state <- workflowSummaryFinalState(config, workflowData, omicType, catFn)
        studyParamsFile <- createWorkflowArgsFn(workflow_name = label, description = description,
            source_dir_path = paths$source_dir, final_s4_object = state$finalS4Object, contrasts_tbl = state$contrastsTbl,
            workflow_data = workflowData,
            omic_type = omicType)
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
            stateName, object, notify = notify, saveRDSFn = saveRDSFn,
            failed = function(error) {
                reasonShown <<- TRUE # the managed step's own notice names the reason
                summaryFailed(error)
            },
            sourceDir = paths$source_dir, inputRemoved = earlier && !dir.exists(inputFolder), done = function(written) {
                # F108-03: a write the gate refused shows no "Saved Integration" notice. KIT-108 W4: a route writer's
                # notice is its own (proteomics': the native archive's, or the S4 RDS's).
                if (!isFALSE(written)) {
                    notify(if (is.character(written)) written else
                        workflowSummaryText(config, "integration_saved", "Saved Integration S4 Object"), "message", 5)
                }
                # KIT-228: no new input (the gate, a writer with no source, or a worker that wrote none): the earlier one
                # was removed at the start; a writer that reported TRUE is judged by its manifest.
                if (isFALSE(written)) notifyInputRemoved(if (isTRUE(attr(written, "gate"))) "gate" else "no_table") else
                    if (!file.exists(file.path(inputFolder %||% "", INTEGRATION_INPUT_MANIFEST))) notifyInputRemoved("no_table")
                notify("Study parameters saved successfully", "message", 5)
                # KIT-216 (SD-2, SD-4): a route writer's reason the session file was not refreshed, after the success.
                if (!is.null(attr(written, "notice_after"))) notify(attr(written, "notice_after"), "warning", NULL)
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
#' @param done `function(written)`: `written` is FALSE when the gate refused the write (with the attribute `gate`, KIT-228:
#'   the Save's notice names the missing export), else TRUE, or a route writer's notice.
#' @param sourceDir The project's source directory, for a route writer (KIT-108 W4: proteomics' native archive).
#' @param inputRemoved Whether this Save removed the omic's earlier integration input at its start (KIT-228).
#' @noRd
startWorkflowSummaryIntegrationExport <- function(config, workflowData, omicType, path, stateName, object, notify, done,
    failed, saveRDSFn = saveRDS, sourceDir = NULL, inputRemoved = FALSE) {
    integration <- config$integration
    requires <- integration$requires
    # KIT-136: a requirement that depends on the session (transcriptomics: the RUVIII export only when RUV was applied).
    if (is.function(requires)) requires <- requires(workflowData)
    if (!is.null(requires) && is.null(workflowData$state_manager$.context$snapshot$exports[[requires]])) {
        notify(integration$missing_message, "warning", 10)
        return(done(structure(FALSE, gate = TRUE)))
    }
    # KIT-082: the integration object is the whole state, so a managed session's worker writes it.
    runner <- if (!is.null(integration$operation)) workflowOmicsManagedActions(omicType)
    if (!is.null(runner) && !is.null(stateName)) {
        # VAL-393: handed to a runner that is retrying a publication, this Save may only finish the earlier request, which
        # does not write the input again: said once.
        if (isTRUE(inputRemoved) && is.function(runner$retrying) && isTRUE(runner$retrying())) {
            notify(workflowSummaryInputRetryText(omicType), "warning", NULL)
        }
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
#' its failures classified by `copy$required` (a required one stops the copy). A combined workbook that could not be
#' read into or saved (A108-04) is reported as `failed` for every omic, whatever `copy$required` names.
#' @param argsSaved Whether the parameters were saved.
#' @return `list(fallbackWritten, required, optional, failed, contrastsTbl, designMatrix, error)`; `error` is the
#'   failure's message, or NULL.
#' @noRd
copyWorkflowSummaryResults <- function(config, workflowData, projectDirs, omicType, label, description, argsSaved, notify,
    copyFn = copyToResultsSummary, catFn = cat) {
    paths <- projectDirs[[omicType]]
    # KIT-212 (VAL-172 (1), (5); E-4): refused first, before the fallback file, when the DA was run on another table or
    # the DA folders hold files no DA of this session made; nothing is written or copied. A catch-all at a press site, as
    # the Report's guard. KIT-243: then refused while the newest DA's publication has not completed (part D), and while the
    # last Save does not describe the current history and DA (VAL-212 (2)).
    refused <- tryCatch({
        assertWorkflowSummaryDaCurrent(workflowData, omicType, paths, files = TRUE)
        assertWorkflowDaPublicationComplete(paths)
        assertWorkflowSummarySaveCurrent(workflowData, paths)
    }, error = function(e) e)
    if (inherits(refused, "error")) {
        notify(paste("Copy error:", cli::ansi_strip(conditionMessage(refused))), "error", 10)
        logger::log_error(logger::skip_formatter(paste("Copy refused:", conditionMessage(refused))))
        return(list(fallbackWritten = FALSE, required = list(), optional = list(), failed = list(), contrastsTbl = NULL,
            designMatrix = NULL, error = conditionMessage(refused)))
    }
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
        # P-30 route B: the copy's workbook transients are freed by R but kept by the allocator, so its free pages are
        # returned once the copy has returned or failed. `copyToResultsSummary` itself is unchanged for notebook callers.
        # KIT-219 (VAL-187): only the files the current enrichment made are copied; NULL (a manager that is not a
        # reference-state session owner) leaves the copy as it was.
        enrichmentFiles <- workflowEnrichmentCurrentFiles(workflowData, omicType)
        # The list changes the copy only when the folder holds a file outside it; with none, the argument is left out and
        # the copy is the unchanged one (R2: the recorded copy calls of the W1 fixtures carry the argument list).
        if (!is.null(enrichmentFiles)) {
            onDisk <- if (isTRUE(dir.exists(paths$pathway_dir %||% ""))) list.files(paths$pathway_dir, recursive = TRUE)
            if (!length(setdiff(onDisk, enrichmentFiles))) enrichmentFiles <- NULL
        }
        failures <- tryCatch(if (is.null(enrichmentFiles)) copyFn(omic_type = omicType, experiment_label = label,
            contrasts_tbl = contrastsTbl, design_matrix = designMatrix, force = TRUE) else copyFn(omic_type = omicType,
            experiment_label = label, contrasts_tbl = contrastsTbl, design_matrix = designMatrix, force = TRUE,
            enrichment_files = enrichmentFiles), finally = reclaimWorkflowSummaryCopyMemory())
        required <- list()
        optional <- list()
        failed <- list()
        notCurrent <- NULL
        if (is.list(failures) && length(failures) > 0L) {
            # The record of the enrichment files left out is not a failure: it is taken out before the classification.
            left <- vapply(failures, function(failure) is.list(failure) && identical(failure$type, "enrichment_not_current"),
                logical(1))
            if (any(left)) notCurrent <- failures[left][[1L]]
            failures <- failures[!left]
        }
        if (is.list(failures) && length(failures) > 0L) {
            workbook <- vapply(failures, workflowSummaryIsWorkbookFailure, logical(1))
            failed <- failures[workbook]
            failures <- failures[!workbook]
        }
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
        # KIT-219: what was left out of the enrichment files, stated after the notice (T1, T1b).
        left <- if (!is.null(notCurrent)) paste0(" ", workflowEnrichmentNotCopiedText(notCurrent, omicType))
        if (length(failed) > 0L) {
            notify(paste0("Publication files copied, but the combined results workbook failed: ",
                workflowSummaryWorkbookFailureText(failed), left), "error", 10)
            logger::log_error(logger::skip_formatter(paste("Combined results workbook failed:",
                workflowSummaryWorkbookFailureText(failed))))
        } else if (!is.null(notCurrent)) {
            notify(paste0(if (length(optional) > 0L) "Publication files copied with optional missing artifacts." else
                "Publication files copied.", left), "warning", 10)
        } else {
            notify(if (length(optional) > 0L) "Publication files copied with optional missing artifacts" else
                "Publication files copied", if (length(optional) > 0L) "warning" else "message", 5)
        }
        list(fallbackWritten = fallbackWritten, required = required, optional = optional, failed = failed,
            contrastsTbl = contrastsTbl, designMatrix = designMatrix, error = NULL)
    }, error = function(e) {
        notify(paste("Copy error:", conditionMessage(e)), "error", 10)
        logger::log_error(logger::skip_formatter(paste("Failed to copy files:", conditionMessage(e))))
        list(fallbackWritten = fallbackWritten, required = list(), optional = list(), failed = list(), contrastsTbl = NULL,
            designMatrix = NULL, error = conditionMessage(e))
    })
}

#' Whether one of the copy's failure records concerns a combined results workbook (A108-04): a DA or enrichment table
#' that could not be read into it, or a workbook that could not be saved (`copyToResultsSummary`'s `da_excel_read`,
#' `enrichment_tsv_read` and `workbook_save` records).
#' @param failure One element of `copyToResultsSummary`'s return value.
#' @return TRUE or FALSE.
#' @noRd
workflowSummaryIsWorkbookFailure <- function(failure) {
    is.list(failure) && isTRUE(failure$type %in% c("workbook_save", "da_excel_read", "enrichment_tsv_read"))
}

#' The workbook failures as one line: each file's name and its error.
#' @param failed The copy's workbook failure records.
#' @return A string.
#' @noRd
workflowSummaryWorkbookFailureText <- function(failed) {
    paste(vapply(failed, function(failure) sprintf("%s (%s)", basename(failure$path %||% failure$source %||% "?"),
        failure$error %||% "failed"), character(1)), collapse = "; ")
}

#' Return the allocator's free pages to the system after the publication copy (P-30 route B): the package's
#' `reclaimWorkflowMemory()` (a full collection, then glibc's `malloc_trim(0)`; the collection alone elsewhere). Optional: a
#' failure is logged and the copy's outcome stands.
#' @return The reclaim's status list, invisibly.
#' @noRd
reclaimWorkflowSummaryCopyMemory <- function() {
    invisible(tryCatch(reclaimWorkflowMemory(), workflow_reclaim_error = function(error) {
        logger::log_warn(logger::skip_formatter(paste("Optional memory reclamation failed:", conditionMessage(error))))
        list(status = "failed")
    }))
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
    # KIT-207 (SD-2 (b)): refused, before the template is looked up, when the head is not the normalisation's result.
    # Recorded divergence: a catch-all at a press site (as KIT-194 A05), over conventions/R.md's specific conditions, so
    # that an unexpected error of the guard's reads is shown rather than escaping the observer.
    # KIT-212 (VAL-172 (1), (5)): then refused when the DA was run on another table or the DA folders hold files no DA of
    # this session made. KIT-243: then while the newest DA's publication has not completed, and while the last Save does
    # not describe the current history and DA.
    refused <- tryCatch({
        assertWorkflowSummaryHeadFinal(workflowData, omicType)
        assertWorkflowSummaryDaCurrent(workflowData, omicType, paths, files = TRUE)
        assertWorkflowDaPublicationComplete(paths)
        assertWorkflowSummarySaveCurrent(workflowData, paths)
        assertWorkflowSummaryEnrichmentCurrent(workflowData, omicType, paths)
    }, error = function(e) e)
    if (inherits(refused, "error")) {
        logger::log_error(logger::skip_formatter(paste("Report refused:", conditionMessage(refused))))
        notify(paste("Report generation failed:", cli::ansi_strip(conditionMessage(refused))), "error", 15)
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
        # KIT-207 (SD-2 (b)): refused when the head is not the normalisation's result ("Export failed: <text>").
        assertWorkflowSummaryHeadFinal(workflowData, omicType)
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
        # KIT-274: the request's setup (the head's reference and the step's parameters, which may refuse) runs inside the
        # failure handler, as the in-process branch's run does.
        actions <- tryCatch({
            source <- getWorkflowStateReference(manager)
            export <- omicsManagedProduct(step$operation, step$parameters(manager, source, choice), function(result) {
                step$accept(manager, source, list(reference = result$reference, value = result$value))
                held$value <- result$value
            })
            c(list(export), lapply(step$publications, function(publication) {
                omicsManagedProduct(publication$operation, publication$parameters, function(result) invisible(NULL))
            }))
        }, error = function(error) {
            failed(error)
            NULL
        })
        if (is.null(actions)) return(invisible(FALSE))
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

#' The integration export in process: KIT-095's writer of the accepted RUVIII matrix (KIT-136: the TMM logCPM when RUV
#' was not applied), then the declared integration input (KIT-136) beside it, at `integration/transcriptomics/`.
#' @noRd
trxSummaryIntegrationWrite <- function(workflowData, path) {
    manager <- workflowData$state_manager
    withWorkflowStateScope(manager, function() trxWriteIntegrationExport(manager, list(), path))
    trxWriteIntegrationInput(transcriptomicsDeclarationBinding, manager, file.path(dirname(path), "transcriptomics"))
}

#' The export the transcriptomics integration needs (the Summary's `integration$requires`, KIT-136): the accepted RUVIII
#' export (INT-1: the corrected matrix), unless the head's RUV step applied no RUV (design `none`, or k of 0), whose
#' input is the TMM logCPM and needs none.
#' @return `"trx_ruviii_export"` or NULL.
#' @noRd
trxSummaryIntegrationRequires <- function(workflowData) {
    manager <- workflowData$state_manager
    if (!inherits(manager, "DIAStateAdapter")) return("trx_ruviii_export")
    state <- withWorkflowStateScope(manager, function() trxIntegrationRuvState(manager))
    if (!is.null(state) && !state$applied) NULL else "trx_ruviii_export"
}

#' The RUVIII export's parameters (D4; `runTrxRuviiiExport`'s at 5be34ee5, mod_trx_summary.R:76-78): the head's RUVs
#' logCPM and grid entities, and M.
#' @param source The head's reference.
#' @noRd
trxRuviiiExportParameters <- function(manager, source, m_source) {
    ruvs <- withWorkflowStateScope(manager, function() {
        # KIT-274: a head that applied no RUV refuses the managed request before anything is sent to the worker.
        trxAssertRuviiiExportApplied(manager)
        trxStateRuvs(manager, source)
    })
    list(logcpm = ruvs$logcpm$entity_id, grid = ruvs$grid$entity_id, ruviii_m_source = m_source)
}

#' The session-state payload: lipid's, with the head's metadata (`trxSummaryStateObject`), the omic and its template
#' (`prepareLipidSummarySessionState` under the transcriptomics profile's entries at 5be34ee5, mod_transcriptomics.R:39,
#' :45-57; KIT-109 passes them as arguments).
#' @noRd
trxSummarySessionState <- function(inputValues, projectDirs, omicType, values, timestamp, workflowData = NULL) {
    prepareLipidSummarySessionState(inputValues = inputValues, projectDirs = projectDirs, omicType = omicType,
        values = values, timestamp = timestamp, workflowData = workflowData, stateObject = trxSummaryStateObject,
        workflowType = "transcriptomics", reportTemplate = "transcriptomics_report.rmd")
}

# ------------------------------------------------------------------------------------------------
# Lipidomics' session-state functions (moved unchanged from mod_lipid_summary_state_helpers.R, D3).
# ------------------------------------------------------------------------------------------------

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
    stateObject = getLipidSummaryStateObject,
    workflowType = "lipidomics",
    reportTemplate = "lipidomics_report.rmd"
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

    finalS4Object <- stateObject(workflowData)
    parameterPayload <- buildLipidSummaryParameterPayload(workflowData, finalS4Object)
    globalParameters <- lipidSummaryFirstNonNull(
        parameterPayload$config_list$globalParameters,
        parameterPayload$s4_args$globalParameters,
        list()
    )

    sessionState$workflow_type <- lipidSummaryFirstNonNull(
        globalParameters$workflow_type,
        workflowType
    )
    sessionState$report_template <- lipidSummaryFirstNonNull(
        globalParameters$report_template,
        reportTemplate
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
