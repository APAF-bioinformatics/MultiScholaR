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

# KIT-138, KIT-139: the integration GUI's configuration, its wording table and the helpers its tabs call (the project
# scan, the import request, the sample-matching preview). The model layer is KIT-062's
# (utils_workflow_integration_model.R); nothing here computes beyond reading what it returns. The whole surface sits
# behind the option `MultiScholaR.integration_gui` (off by default): with it off, the home page and the app are today's,
# byte for byte. The single-omic `validateOmicGuiConfig` does not fit (it requires design, qc and da groups); the shell
# is validated by `validateOmicGuiShell` only (KIT-060 audit section 5.3).

#' Whether the integration GUI is offered (KIT-138; the gate in the manner of `MultiScholaR.ruv_design_panel`).
#' `options(MultiScholaR.integration_gui = TRUE)` shows it.
#' @return One logical.
#' @noRd
integrationGuiEnabled <- function() {
    isTRUE(getOption("MultiScholaR.integration_gui", FALSE))
}

#' Every user-facing string of the integration GUI (KIT-138, KIT-139): one table, drafts the user sets the final
#' wording of (docs/qa/workflow-execution/kit/kit-138/wording.md). The RUV approaches' names are the RUV panel's
#' (`omicGuiRuvRecordTexts`), so that one approach reads the same everywhere. `%s`/`%d` are filled in the order the
#' comment beside each entry gives.
#' @return A named character vector.
#' @noRd
integrationGuiTexts <- function() {
    ruv <- omicGuiRuvRecordTexts()
    c(
        # The home page (gate on).
        home_tile_title = "Integration",
        home_tile_text = "Combine two or more omics already processed in this project",
        # The tabs and the stepper.
        tab_import = " Import",
        step_import = "Import",
        tab_alignment = " Sample Matching",
        step_alignment = "Samples",
        workflow_label = "Integration",
        # The Import tab.
        import_heading = "Import processed omics",
        import_intro = paste("MultiScholaR looks in this project for the omics saved for integration at their Summary",
            "tab. Choose the data sets to combine; two or more omics are needed."),
        # %s: the folder.
        import_folder = "Folder searched: %s",
        import_rescan = "Search again",
        # %d: ready, then not ready.
        import_found = "Omics ready to combine: %d. Not ready: %d.",
        import_none = paste("No omic saved for integration was found in this folder. Open each omic's Summary tab, save",
            "its results for integration, then search again."),
        omics_heading = "Omics found",
        col_omic = "Omic",
        col_data_type = "Data type",
        col_values = "Values",
        col_correction = "Correction for unwanted variation (RUV)",
        col_saved = "Saved",
        col_status = "Status",
        values_corrected = "corrected for unwanted variation",
        values_normalised = "normalised, not corrected",
        correction_none = "not applied",
        # %s: the approach.
        correction_design = "applied: %s",
        # %s: the method, then the replicate sets.
        correction_legacy = "applied: the approach of earlier versions (%s; replicate sets: %s)",
        # %s: the number of factors.
        correction_k = "%s factor(s)",
        # %s: the data set, then the number of factors.
        correction_view_applied = "%s: %s factor(s)",
        # %s: the data set.
        correction_view_not_applied = "%s: not applied",
        correction_unknown = "applied (approach not recorded)",
        method_ruviii_c = "RUV-III-C",
        method_ruviii = "RUV-III",
        sets_group = ruv[["sets_group"]],
        sets_technical_replicates = ruv[["sets_technical_replicates"]],
        design_groups_all_features = ruv[["design_groups_all_features"]],
        design_technical_replicates = ruv[["design_technical_replicates"]],
        design_supplied_controls = ruv[["design_supplied_controls"]],
        design_first_pass_nulls = ruv[["design_first_pass_nulls"]],
        design_none = ruv[["design_none"]],
        # %s: the time of saving.
        saved_at = "%s UTC",
        status_ready = "ready",
        status_changed = "changed since it was imported: import again",
        status_not_saved = "not saved for integration: save it at that omic's Summary tab",
        status_refused = "cannot be used",
        refused_heading = "Not ready",
        # %s: the omic, then the reason (the check's own message).
        refused_reason = "%s: %s",
        views_heading = "Data sets found",
        col_view = "Data set",
        col_features = "Features",
        col_runs = "Runs",
        col_biological_samples = "Biological samples",
        col_groups = "Groups",
        view_choice_label = "Data sets to combine",
        # %s: the data set, the omic; %d: features, biological samples.
        view_choice_item = "%s (%s; %d features, %d biological samples)",
        import_button = "Import selected data",
        # %d: the omics selected.
        import_need_two_omics = "Select data sets from two or more omics; %d selected.",
        # %d: the data sets selected.
        import_need_two_views = "Select two or more data sets; %d selected.",
        import_need_ready = "Only omics marked ready can be imported.",
        # %s: the reason.
        import_refused = "The import was refused: %s",
        # %d: omics, data sets, runs, biological samples.
        import_done = "Imported %d omics: %d data sets, %d runs, %d biological samples.",
        # %s: the reason.
        store_unavailable = "Integration cannot keep its work in this project: %s",
        # The Sample Matching tab.
        align_heading = "Match samples across omics",
        align_intro = paste("Each row is one biological sample, matched across omics by the sample column each omic",
            "saved. Check that every row is one biological sample, then confirm the table."),
        align_not_imported = "Import two or more omics at the Import tab first.",
        align_key_heading = "How samples are matched",
        # %s: the omic, then its column.
        align_key = "%s: by the column %s",
        align_caution = paste("Samples are matched by name: a sample must carry the same name in every omic. A sample",
            "named differently in two omics appears as two rows."),
        align_rules_heading = "What happens before the model",
        align_provisional = "These rules are provisional: the integration methods may still change them.",
        # %d: all biological samples, those measured in every data set.
        fact_samples = "%d biological samples in all; %d measured in every data set.",
        # %d: samples missing from at least one data set.
        fact_missing_keep_missing = paste("%d biological samples are missing from at least one data set; they are kept,",
            "with that data set's values missing."),
        fact_missing_drop_incomplete = "%d biological samples missing from at least one data set are left out.",
        # %s: the data set; %d: runs, biological samples, samples measured more than once.
        fact_replicates_average = paste("%s: %d runs from %d biological samples; the runs of the %d samples measured more",
            "than once are averaged."),
        # %s: the data set; %d: runs. Every biological sample measured once.
        fact_replicates_single = "%s: %d runs, one per biological sample; nothing to average.",
        # %s: the data set; %d: internal standards.
        fact_standards_drop = "%s: %d internal standards removed.",
        fact_standards_keep = "%s: %d internal standards kept.",
        # %s: the data set; %d: features.
        fact_sparse = "%s: %d features measured in fewer than two samples are left out.",
        # %d: samples in every data set, the threshold.
        few_samples = paste("Only %d biological samples are measured in every data set (fewer than %d). A factor model",
            "needs more samples; its results may be unstable."),
        view_table_heading = "Data sets",
        col_standards_removed = "Internal standards removed",
        col_sparse = "Features left out (fewer than two samples)",
        col_features_used = "Features used",
        sample_table_heading = "Biological samples",
        col_sample = "Biological sample",
        col_group = "Group",
        col_views_present = "Data sets measured",
        cell_missing = "missing",
        confirm_button = "Confirm this sample table",
        # %s: the time, then the file.
        confirmed = "Sample table confirmed at %s and saved to %s.",
        not_confirmed = "Not confirmed yet.",
        # %s: the reason.
        align_refused = "The samples cannot be matched: %s",
        # %s: the data sets.
        refuse_no_shared = "These data sets share no biological sample (%s). Check that each omic names its samples the same way.",
        # %d: the data sets selected.
        refuse_too_few_views = "Two or more data sets are needed; %d selected.",
        # KIT-140: the check's own message is shown under "Technical detail", not in the sentence.
        refuse_changed = "A data set has changed since it was imported. Import the data again at the Import tab.",
        # KIT-140: the refusals in plain words, by the refusal's condition class (`integrationRefusalText`); the check's
        # own message is kept under a collapsed "Technical detail" line.
        technical_detail = "Technical detail",
        refusal_input_invalid = paste("This omic's saved files cannot be used as they are: they are incomplete, or were",
            "changed after they were saved. Save the omic again at its Summary tab, then search again."),
        refusal_input_stale = "A data set has changed since it was imported. Import the data again at the Import tab.",
        refusal_input_duplicate = paste("Two inputs are of the same omic, or two data sets have the same name. Keep one",
            "input per omic."),
        refusal_too_few_omics = "Two or more omics are needed.",
        refusal_too_few_views = "Two or more data sets are needed.",
        refusal_inputs_missing = "Nothing is imported yet. Import two or more omics at the Import tab first.",
        refusal_no_shared = paste("The data sets share no biological sample. Check that each omic names its samples the",
            "same way."),
        refusal_group_conflict = paste("A biological sample is in different groups in different omics. Check the groups",
            "in each omic's design."),
        refusal_view_empty = "A data set has no feature measured in two or more samples, so nothing of it is left for the model.",
        refusal_parameters_invalid = "A model setting is not valid. Check the values entered.",
        refusal_mofa_unavailable = paste("The model cannot be trained on this computer: the MOFA2 software, or the Python",
            "it uses, is not installed."),
        refusal_other = "This step could not be completed; the technical detail below says why."
    )
}

#' The integration GUI's configuration (KIT-138): the shell group only, in the single-omic shell's form, validated by
#' `validateOmicGuiShell`. KIT-140 adds its tab to `tabs`.
#' @return `list(shell, texts)`.
#' @noRd
integrationGuiConfig <- function() {
    texts <- integrationGuiTexts()
    tab <- function(status_key, value, title, icon, step) {
        omicGuiTab(status_key, value, value, paste0("integration_", value), texts[[title]], icon,
            paste0("integration-tab-", value), texts[[step]])
    }
    tabs <- list(
        tab("integration_import", "import", "tab_import", "file-import", "step_import"),
        tab("integration_alignment", "alignment", "tab_alignment", "th", "step_alignment"))
    shell <- list(omic_type = "integration", prefix = "integration", testid = "integration",
        label = texts[["workflow_label"]], label_lower = "integration", header_icon = NULL,
        tabset_id = "integration_tabs", tab_style = "icon_in_title", routes = "integration", tabs = tabs,
        bus = list(state_manager = NULL, tab_status = NULL, integration_scan = NULL, integration_import = NULL,
            integration_alignment = NULL),
        bus_registration = "userData")
    validateOmicGuiShell(shell, omicGuiConfigSchema())
    list(shell = shell, texts = texts)
}

#' The home page's Integration tile with the gate on: selectable like an omic (the omics it combines were processed
#' earlier, so no other omic needs to be selected in this session).
#' @noRd
integrationHomeTile <- function(texts = integrationGuiTexts()) {
    testid(
        shiny::div(
            id = "integration_box",
            class = "omic-selection-box",
            onclick = "toggleOmicSelection('integration')",
            shiny::h3(shiny::icon("project-diagram"), texts[["home_tile_title"]]),
            shiny::p(texts[["home_tile_text"]])
        ),
        "tile-integration"
    )
}

#' The home page's selection script: as it is with the gate off; with the gate on, Integration is selected like an
#' omic and its tile is never disabled. Refused when the script no longer holds the lines it changes.
#' @param script The home page's script.
#' @return The script.
#' @noRd
integrationHomeScript <- function(script) {
    if (!integrationGuiEnabled()) return(script)
    edits <- c("if (omic === 'integration') return; // Integration is auto-enabled" =
            "// Integration is selected like an omic (KIT-138)",
        "if (selectedOmics.length >= 2) {" = "if (true) { // KIT-138: the Integration tile stays enabled")
    for (from in names(edits)) {
        if (!grepl(from, script, fixed = TRUE)) {
            cli::cli_abort("the home page script no longer holds {.val {from}}", class = "integration_gui_invalid")
        }
        script <- sub(from, edits[[from]], script, fixed = TRUE)
    }
    script
}

#' The "module not available" panel app_server shows for an omic without a module (the integration module with the
#' gate off shows the same panel, so the app is as before).
#' @noRd
integrationGuiUnavailablePanel <- function() {
    shiny::div(
        class = "alert alert-warning",
        shiny::h4("Module not yet implemented"),
        shiny::p(paste("The", "integration", "workflow module is under development."))
    )
}

# ------------------------------------------------------------------------------------------------
# The Import tab: the project scan (KIT-138).
# ------------------------------------------------------------------------------------------------

INTEGRATION_GUI_OMICS <- c("lipidomics", "metabolomics", "proteomics", "transcriptomics")

#' One value of a manifest field read by `jsonlite::read_json(simplifyVector = FALSE)`: NA when absent, null or empty.
#' @noRd
integrationGuiScalar <- function(value) {
    if (is.null(value) || !length(value)) return(NA)
    value <- unlist(value)
    if (!length(value)) NA else value[[1L]]
}

#' The correction for unwanted variation an input names, in words (KIT-138): not applied for a normalised matrix; for
#' a corrected one the declared RUV approach (KIT-129's `ruv_design`), else the approach of earlier versions (its method
#' and replicate sets), then the number of factors, overall or per data set (lipidomics, metabolomics).
#' @param manifest The input's manifest.
#' @param texts `integrationGuiTexts()`.
#' @return One string.
#' @noRd
integrationCorrectionText <- function(manifest, texts = integrationGuiTexts()) {
    correction <- manifest$correction
    if (!identical(manifest$matrix_role, "corrected") || !is.list(correction)) return(texts[["correction_none"]])
    design <- integrationGuiScalar(correction$ruv_design)
    method <- integrationGuiScalar(correction$method)
    sets <- integrationGuiScalar(correction$replicate_structure)
    key <- paste0("design_", design)
    text <- if (!is.na(design) && key %in% names(texts)) {
        sprintf(texts[["correction_design"]], texts[[key]])
    } else if (!is.na(method)) {
        method_text <- switch(as.character(method), RUVIII_C_Varying = texts[["method_ruviii_c"]],
            RUVIII = texts[["method_ruviii"]], as.character(method))
        sets_text <- if (identical(sets, "technical_replicates")) texts[["sets_technical_replicates"]] else texts[["sets_group"]]
        sprintf(texts[["correction_legacy"]], method_text, sets_text)
    } else {
        texts[["correction_unknown"]]
    }
    if (is.list(correction$views) && length(correction$views)) {
        per_view <- vapply(correction$views, function(view) {
            assay <- as.character(integrationGuiScalar(view$assay))
            if (isTRUE(as.logical(integrationGuiScalar(view$applied)))) {
                sprintf(texts[["correction_view_applied"]], assay, as.character(integrationGuiScalar(view$k)))
            } else {
                sprintf(texts[["correction_view_not_applied"]], assay)
            }
        }, character(1))
        return(paste0(text, "; ", paste(per_view, collapse = "; ")))
    }
    k <- integrationGuiScalar(correction$ruv_k)
    if (is.na(k)) k <- integrationGuiScalar(correction$k)
    if (!is.na(k)) text <- paste0(text, "; ", sprintf(texts[["correction_k"]], as.character(k)))
    text
}

#' Scan a project's integration folder for the omics' integration inputs (KIT-138). Each omic's folder
#' (`<project>/integration/<omic>/`) is read by the contract's reader (`readIntegrationInput`): a folder without the
#' input is not saved; one the reader refuses is listed with the reader's message; one whose manifest has changed
#' since the import is marked changed. Nothing is written.
#' @param integration_dir The project's integration folder.
#' @param imported The imported manifests' sha256, named by omic (NULL before an import).
#' @param texts `integrationGuiTexts()`.
#' @return `list(dir, omics, views)`: one row per omic folder found (`omic`, `dir`, `status`: `ready`, `changed`,
#'   `not_saved` or `refused`; `reason`, `data_type`, `values`, `correction`, `saved`, `manifest_sha256`) and one row
#'   per data set of a usable omic (`view`, `omic`, `features`, `runs`, `biological_samples`, `groups`).
#' @noRd
integrationScanProject <- function(integration_dir, imported = NULL, texts = integrationGuiTexts()) {
    folders <- file.path(integration_dir, INTEGRATION_GUI_OMICS)
    folders <- folders[dir.exists(folders)]
    omics <- list()
    views <- list()
    for (dir in folders) {
        omic <- basename(dir)
        row <- list(omic = omic, dir = normalizePath(dir, winslash = "/", mustWork = TRUE), status = "ready",
            reason = NA_character_, data_type = NA_character_, values = NA_character_, correction = NA_character_,
            saved = NA_character_, manifest_sha256 = NA_character_)
        manifest_path <- file.path(dir, INTEGRATION_INPUT_MANIFEST)
        if (!file.exists(manifest_path)) {
            row$status <- "not_saved"
            row$reason <- texts[["status_not_saved"]]
            omics[[omic]] <- row
            next
        }
        input <- tryCatch(readIntegrationInput(dir), integration_input_invalid = function(error) error)
        if (inherits(input, "error")) {
            row$status <- "refused"
            row$reason <- conditionMessage(input)
            omics[[omic]] <- row
            next
        }
        manifest <- input$manifest
        row$manifest_sha256 <- digest::digest(file = manifest_path, algo = "sha256")
        if (!is.null(imported) && omic %in% names(imported) && !identical(unname(imported[[omic]]), row$manifest_sha256)) {
            row$status <- "changed"
            row$reason <- texts[["status_changed"]]
        }
        row$data_type <- as.character(integrationGuiScalar(manifest$workflow))
        row$values <- texts[[if (identical(manifest$matrix_role, "corrected")) "values_corrected" else "values_normalised"]]
        row$correction <- integrationCorrectionText(manifest, texts)
        written <- as.character(integrationGuiScalar(manifest$written_at))
        row$saved <- if (is.na(written)) NA_character_ else
            sprintf(texts[["saved_at"]], sub("^(\\d{4}-\\d{2}-\\d{2})T(\\d{2}:\\d{2}).*$", "\\1 \\2", written))
        omics[[omic]] <- row
        samples <- input$samples
        biological <- manifest$biological_sample_column %||% integrationInputBiologicalSampleColumn(samples,
            integrationGuiScalar(manifest$technical_replicate_column), manifest$sample_id_column)
        biological_count <- if (is.character(biological) && biological %in% names(samples))
            length(unique(samples[[biological]])) else NA_integer_
        group_count <- length(unique(samples[[manifest$group_column]]))
        for (view in names(input$views)) {
            matrix <- input$views[[view]]$matrix
            views[[view]] <- data.frame(view = view, omic = omic, features = nrow(matrix), runs = ncol(matrix) - 1L,
                biological_samples = biological_count, groups = group_count, stringsAsFactors = FALSE)
        }
    }
    omic_table <- if (length(omics)) do.call(rbind, lapply(omics, as.data.frame, stringsAsFactors = FALSE)) else
        data.frame(omic = character(), dir = character(), status = character(), reason = character(),
            data_type = character(), values = character(), correction = character(), saved = character(),
            manifest_sha256 = character(), stringsAsFactors = FALSE)
    rownames(omic_table) <- NULL
    view_table <- if (length(views)) do.call(rbind, views) else data.frame(view = character(), omic = character(),
        features = integer(), runs = integer(), biological_samples = integer(), groups = integer(),
        stringsAsFactors = FALSE)
    rownames(view_table) <- NULL
    list(dir = integration_dir, omics = omic_table, views = view_table)
}

#' The scan's omics as the Import tab shows them: one row per omic folder, the headers in words.
#' @noRd
integrationScanDisplay <- function(scan, texts = integrationGuiTexts()) {
    omics <- scan$omics
    status <- vapply(seq_len(nrow(omics)), function(i) switch(omics$status[[i]], ready = texts[["status_ready"]],
        changed = texts[["status_changed"]], not_saved = texts[["status_not_saved"]], texts[["status_refused"]]),
        character(1))
    table <- data.frame(tools::toTitleCase(omics$omic), omics$data_type, omics$values, omics$correction, omics$saved,
        status, stringsAsFactors = FALSE)
    names(table) <- unname(texts[c("col_omic", "col_data_type", "col_values", "col_correction", "col_saved", "col_status")])
    table
}

#' The data sets of the omics an import can take: those of the omics marked ready, or changed since the import (a new
#' import is the remedy).
#' @noRd
integrationScanReadyViews <- function(scan) {
    ready <- scan$omics$omic[scan$omics$status %in% c("ready", "changed")]
    scan$views[scan$views$omic %in% ready, , drop = FALSE]
}

#' The import a selection asks for (KIT-138): the folders of the omics with a selected data set and the selected data
#' sets in the scan's order. Refused in words (`integration_refused`): a data set of an omic not ready, fewer than two
#' omics (D7), fewer than two data sets.
#' @param scan `integrationScanProject()`'s value.
#' @param selected The selected data sets.
#' @return `list(dirs, omics, views)`.
#' @noRd
integrationImportRequest <- function(scan, selected, texts = integrationGuiTexts()) {
    selected <- unique(as.character(selected))
    ready <- integrationScanReadyViews(scan)
    refuse <- function(message, class) integrationRefuse(message, c(class, "integration_gui_refusal"))
    if (!all(selected %in% ready$view)) refuse(texts[["import_need_ready"]], "integration_input_not_ready")
    views <- ready$view[ready$view %in% selected]
    omics <- unique(ready$omic[ready$view %in% views])
    if (length(omics) < 2L) refuse(sprintf(texts[["import_need_two_omics"]], length(omics)), "integration_too_few_omics")
    if (length(views) < 2L) refuse(sprintf(texts[["import_need_two_views"]], length(views)), "integration_too_few_views")
    list(dirs = scan$omics$dir[match(omics, scan$omics$omic)], omics = omics, views = views)
}

#' Commit an import on the integration store (KIT-138): `importIntegrationInputs` on the folders, then the summary the
#' tabs read (the state's reference, the omics, the selected data sets in the state's order, the manifests' sha256 and
#' the counts).
#' @param manager The integration store (`newIntegrationStateOwner`).
#' @param request `integrationImportRequest()`'s value.
#' @return The import's summary.
#' @noRd
integrationCommitImport <- function(manager, request) {
    reference <- importIntegrationInputs(manager, request$dirs)
    object <- integrationImportedObject(manager)
    views <- names(object@views)[names(object@views) %in% request$views]
    design <- object@design_matrix
    view_omics <- unlist(unname(lapply(object@inputs, function(entry) stats::setNames(rep(entry$omic, length(entry$views)),
        entry$views))))
    used <- design$omic %in% unique(view_omics[views])
    list(reference = reference, omics = names(object@inputs), views = views, manifests = integrationStateManifests(object),
        runs = sum(used), biological_samples = length(unique(design$biological_sample[used])), committed_at = Sys.time())
}

#' The imported inputs state, every field hydrated.
#' @noRd
integrationImportedObject <- function(manager) {
    withWorkflowStateScope(manager, function() {
        source <- getWorkflowStateReference(manager, INTEGRATION_INPUTS_STATE)
        readWorkflowExecutionComponents(manager, source, fields = TRUE)
    })
}

# ------------------------------------------------------------------------------------------------
# The Sample Matching tab (KIT-139).
# ------------------------------------------------------------------------------------------------

#' The sample-matching preview of the imported data sets (KIT-139): the model layer's alignment
#' (`integrationAlignViews`) under the parameters as resolved for the selected data sets (the provisional defaults of
#' `integrationMofaParameterSpec`, never set here), with its facts counted for the tab. The few-samples warning (D7)
#' is kept as a note, not raised. Refused (`integration_refused`): an input changed since the import, fewer than two
#' data sets, no shared biological sample, a sample in two groups.
#' @param object The imported inputs state (`integrationImportedObject`).
#' @param views The selected data sets.
#' @return `list(parameters, aligned, few_samples, keys, views, replicates)`.
#' @noRd
integrationAlignmentPreview <- function(object, views) {
    integrationAssertInputsCurrent(object)
    parameters <- integrationResolveMofaParameters(object, list(views = views))
    few <- FALSE
    aligned <- withCallingHandlers(integrationAlignViews(object, parameters), integration_few_samples = function(warning) {
        few <<- TRUE
        invokeRestart("muffleWarning")
    })
    design <- object@design_matrix
    view_omic <- unlist(unname(lapply(object@inputs, function(entry) stats::setNames(rep(entry$omic, length(entry$views)),
        entry$views))))
    all_samples <- character()
    replicates <- do.call(rbind, lapply(parameters$views, function(view) {
        runs <- names(object@views[[view]])[-1L]
        rows <- design[design$omic == view_omic[[view]] & design$run %in% runs, , drop = FALSE]
        all_samples <<- union(all_samples, rows$biological_sample)
        per_sample <- table(rows$biological_sample)
        standards <- sum(integrationInternalStandards(object@features[[view]]))
        data.frame(view = view, omic = view_omic[[view]], runs = length(runs), biological_samples = length(per_sample),
            repeated = sum(per_sample > 1L), standards = standards, stringsAsFactors = FALSE)
    }))
    keys <- data.frame(omic = names(object@inputs), column = vapply(object@inputs, function(entry)
        as.character(entry$biological_sample_column), character(1)), stringsAsFactors = FALSE, row.names = NULL)
    list(parameters = parameters, aligned = aligned, few_samples = few, keys = keys, views = parameters$views,
        replicates = replicates, all_samples = length(all_samples))
}

#' The facts the Sample Matching tab states, in words, each from the preview's counts and the parameter's value (so a
#' changed default changes the sentence, not a hard-coded rule).
#' @param preview `integrationAlignmentPreview()`'s value.
#' @return A character vector.
#' @noRd
integrationAlignmentFacts <- function(preview, texts = integrationGuiTexts()) {
    parameters <- preview$parameters
    summary <- preview$aligned$view_summary
    text <- function(prefix, value) texts[[paste0(prefix, value)]]
    # Counted over the union of the data sets' biological samples, whichever rule keeps them.
    total <- preview$all_samples
    missing <- total - preview$aligned$shared_samples
    facts <- c(sprintf(texts[["fact_samples"]], total, preview$aligned$shared_samples),
        sprintf(text("fact_missing_", parameters$missing_samples), missing))
    for (i in seq_len(nrow(preview$replicates))) {
        row <- preview$replicates[i, ]
        facts <- c(facts, if (row$repeated == 0L) sprintf(texts[["fact_replicates_single"]], row$view, row$runs) else
            sprintf(text("fact_replicates_", parameters$technical_replicates), row$view, row$runs,
                row$biological_samples, row$repeated))
        if (row$standards > 0L) {
            facts <- c(facts, sprintf(text("fact_standards_", parameters$internal_standards), row$view, row$standards))
        }
        sparse <- summary$sparse_dropped[summary$view == row$view]
        if (length(sparse) && sparse > 0L) facts <- c(facts, sprintf(texts[["fact_sparse"]], row$view, sparse))
    }
    unname(facts)
}

#' The Sample Matching tab's table of data sets: runs, biological samples, internal standards removed, features left
#' out and used.
#' @noRd
integrationAlignmentViewTable <- function(preview, texts = integrationGuiTexts()) {
    summary <- preview$aligned$view_summary
    replicates <- preview$replicates[match(summary$view, preview$replicates$view), , drop = FALSE]
    table <- data.frame(summary$view, replicates$runs, replicates$biological_samples, summary$internal_standards_dropped,
        summary$sparse_dropped, summary$features_model, stringsAsFactors = FALSE)
    names(table) <- unname(texts[c("col_view", "col_runs", "col_biological_samples", "col_standards_removed",
        "col_sparse", "col_features_used")])
    table
}

#' The Sample Matching tab's sample table: one row per biological sample, its group, the number of data sets it is
#' measured in, and per data set the runs it comes from or "missing" (D4).
#' @noRd
integrationAlignmentSampleTable <- function(preview, texts = integrationGuiTexts()) {
    samples <- preview$aligned$samples
    table <- data.frame(samples$biological_sample, samples$group, samples$views_present, stringsAsFactors = FALSE)
    names(table) <- unname(texts[c("col_sample", "col_group", "col_views_present")])
    for (view in preview$views) {
        runs <- samples[[paste0(view, "_runs")]]
        table[[view]] <- ifelse(nzchar(runs), gsub(";", ", ", runs, fixed = TRUE), texts[["cell_missing"]])
    }
    table
}

#' A refusal of the sample matching in words: the tab's own sentence for the refusals it names (KIT-139: no shared
#' sample, fewer than two data sets, a changed input), the plain sentence of the refusal's class for the others
#' (KIT-140, `integrationRefusalText`). The check's own message is the refusal's detail, shown apart.
#' @param condition The `integration_refused` condition.
#' @param views The selected data sets.
#' @noRd
integrationAlignmentRefusal <- function(condition, views, texts = integrationGuiTexts()) {
    reason <- if (inherits(condition, "integration_no_shared_samples")) {
        sprintf(texts[["refuse_no_shared"]], paste(views, collapse = ", "))
    } else if (inherits(condition, "integration_too_few_views")) {
        sprintf(texts[["refuse_too_few_views"]], length(views))
    } else if (inherits(condition, "integration_input_stale")) {
        texts[["refuse_changed"]]
    } else {
        integrationRefusalText(condition, texts)$sentence
    }
    sprintf(texts[["align_refused"]], reason)
}

# ------------------------------------------------------------------------------------------------
# Refusals in plain words (KIT-140): the GUI's sentence by the refusal's condition class. The classes and the messages
# the checks raise are unchanged (tests and records rely on them); the message is kept as the refusal's detail.
# ------------------------------------------------------------------------------------------------

#' The plain sentence of each refusal class the integration tabs can meet, in the order they are looked up (a refusal
#' carries its own class and `integration_refused`).
#' @noRd
INTEGRATION_GUI_REFUSALS <- c(
    integration_input_invalid = "refusal_input_invalid",
    integration_input_stale = "refusal_input_stale",
    integration_input_duplicate = "refusal_input_duplicate",
    integration_too_few_omics = "refusal_too_few_omics",
    integration_too_few_views = "refusal_too_few_views",
    integration_inputs_missing = "refusal_inputs_missing",
    integration_no_shared_samples = "refusal_no_shared",
    integration_group_conflict = "refusal_group_conflict",
    integration_view_empty = "refusal_view_empty",
    integration_parameters_invalid = "refusal_parameters_invalid",
    integration_mofa_unavailable = "refusal_mofa_unavailable")

#' A refusal in words: a sentence the tab wrote itself (`integration_gui_refusal`) as it is; any other refusal as the
#' plain sentence of its class, the check's own message as its detail.
#' @param condition The condition.
#' @param texts `integrationGuiTexts()`.
#' @return `list(sentence, detail)`; `detail` NULL for the tab's own sentences.
#' @noRd
integrationRefusalText <- function(condition, texts = integrationGuiTexts()) {
    if (inherits(condition, "integration_gui_refusal")) return(list(sentence = conditionMessage(condition), detail = NULL))
    known <- names(INTEGRATION_GUI_REFUSALS)[vapply(names(INTEGRATION_GUI_REFUSALS), function(class) inherits(condition, class),
        logical(1))]
    key <- if (length(known)) INTEGRATION_GUI_REFUSALS[[known[[1L]]]] else "refusal_other"
    list(sentence = texts[[key]], detail = conditionMessage(condition))
}

#' Refuse with one of the tab's own sentences (KIT-140): a plain condition, so that the sentence reaches the screen
#' exactly as written (no console formatting); classed `integration_gui_refusal` and `integration_refused`.
#' @param sentence The sentence.
#' @param class The refusal's own class.
#' @noRd
integrationGuiRefuse <- function(sentence, class) {
    stop(structure(class = c(class, "integration_gui_refusal", "integration_refused", "error", "condition"),
        list(message = sentence, call = NULL)))
}

#' A refusal on screen: the sentence, and the check's message under a collapsed "Technical detail" line.
#' @param sentence The sentence.
#' @param detail The check's message, or NULL.
#' @param class The box's class.
#' @noRd
integrationRefusalUi <- function(sentence, detail = NULL, class = "alert alert-danger", texts = integrationGuiTexts()) {
    shiny::div(class = class, shiny::p(sentence),
        if (!is.null(detail) && length(detail) && nzchar(detail)) {
            shiny::tags$details(shiny::tags$summary(texts[["technical_detail"]]), shiny::tags$code(detail))
        })
}

#' Save the confirmed sample table as the integration design (KIT-139): the aligned sample table, tab-separated, in the
#' integration results' inputs folder.
#' @param samples The aligned sample table (`integrationAlignViews()$samples`).
#' @param dir The folder.
#' @return The file's path.
#' @noRd
integrationSaveDesign <- function(samples, dir) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    path <- file.path(dir, "integration_design.tsv")
    vroom::vroom_write(samples, path, delim = "\t", na = "NA", progress = FALSE)
    path
}
